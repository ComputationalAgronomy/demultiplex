import gzip
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from io import BufferedWriter
import logging
import multiprocessing as mp
import os
from rapidfuzz import process
from rapidfuzz.distance import Levenshtein
import shutil
from tqdm import tqdm


FORMAT = "%(message)s"
logger = logging.getLogger()
logger.setLevel(logging.INFO)
formatter = logging.Formatter(FORMAT)

sh = logging.StreamHandler()
sh.setLevel(logging.INFO)
sh.setFormatter(formatter)
logger.addHandler(sh)


class FuzzyMatch():
    def __init__(self, choices: dict[str, str]):
        self.choices = choices

    def best_levenshtein(self, query, allow_dist):
        return process.extractOne(query, self.choices, scorer=Levenshtein.distance, score_cutoff=allow_dist)

    def all_levenshtein(self, query, allow_dist):
        return process.extract(query, self.choices, scorer=Levenshtein.distance, score_cutoff=allow_dist)


class Demux():
    def __init__(self, save_dir, allow_mismatch, threads):
        self.save_dir = save_dir
        self.allow_mismatch = allow_mismatch
        self._get_threads(threads)

    def _get_threads(self, threads):
        max = mp.cpu_count() - 1
        if threads is None:
            logger.info(f"threads wasn't specified. Using maximum available threads {max}.")
            self.threads = max
        else:
            try:
                threads = int(threads)
            except ValueError:
                logger.warning(f"Invalid value for threads. Using maximum available threads {max}.")
                self.threads = max

            if 0 < threads <= max:
                self.threads = threads
            else:
                logger.warning(f"Threads must be between 1 and {max}. Using maximum available threads {max}.")
                self.threads = max

    @staticmethod
    def _approx_match(seq, fm, allow_dist) -> str:
        result = fm.best_levenshtein(seq, allow_dist)
        if result:
            return result[2]
        return None

    def _get_chunks(self, fq_gz):
        records = []
        with gzip.open(fq_gz, 'rt') as in_handle:
            for title, seq, _ in tqdm(FastqGeneralIterator(in_handle), desc="Reading raw data..."):
                seq_id = title.split(" ")[0]
                records.append((seq_id, seq))

        chunk_size = len(records) // self.threads + 1
        self.chunks = [records[i:i+chunk_size] for i in range(0, len(records), chunk_size)]

    @staticmethod
    def _rm_empty(match_dict):
        return {k: v for k, v in match_dict.items() if v!=[]}

    @staticmethod
    def _get_records(fq_gz, fq):
        with gzip.open(fq_gz, 'rb') as f_in:
            with open(fq, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        records = SeqIO.index(fq, "fastq")
        return records

    @staticmethod
    def _write_fq(records, query_list, save_path):
        with BufferedWriter(gzip.open(save_path, 'ab')) as out_handle:
            for seq_id in tqdm(query_list, total=len(query_list), desc=f"Writing {save_path}"):
                out_handle.write(records.get_raw(seq_id))

    @staticmethod
    def _del_fq(records, fa):
        del records
        os.remove(fa)


class SingleEndDemux(Demux):
    def __init__(self,
            index_path,
            raw_path,
            save_dir,
            allow_mismatch: int = 0,
            threads: int = None,
            read_direction: str = "forward",
            pr: str = "GTCGGTAAAACTCGTGCCAGC",
            dry: bool = False
        ):
        """
        Demultiplex single-end fastq.gz file to different sample fastq.gz files.

        :param index_path: path to barcode file, format: <sample_id>\\t<barcode>.
        :param raw_path: path to raw .fastq.gz file.
        :param save_dir: The directory to save demultiplexed fastq.gz file.
        :param pr: The primer sequence after index. Default is "GTCGGTAAAACTCGTGCCAGC"(MiFish-UF).
        :param allow_mismatch: allow distance for barcode matching. If not specify, automatically use proper number. Default is None
        :param threads: Number of CPU cores to be used for processing. If not specify, use all. Default is None.
        :param read_direction: forward or reverse. Default is "forward"
        :param dry: Dry run. If True, no files will be written; only a demultiplexed report will be output. Default is False.
        """
        super().__init__(save_dir, allow_mismatch, threads)
        self.suffix = "R1" if read_direction in ["forward", "Forward", "F", "R1"] else "R2"
        self.match_dict = {}
        self.undetermined = []
        self._read_barcode(index_path)
        self._init_dict()
        self.fm = FuzzyMatch(self.barcode)
        self.pr = pr
        self.run(raw_path, dry)

    def _read_barcode(self, barcode_path):
        try:
            self.sample_id_list = []
            self.barcode = {}
            with open(barcode_path, 'r') as f:
                for line in tqdm(f, desc="Reading Index file..."):
                    line = line.strip().split("\t")
                    self.sample_id_list.append(line[0])
                    self.barcode[line[0]] = line[1]
            self.barcode_len = len(line[1])
        except Exception as e:
            logger.error(f"Error reading Index file: {e}. Error line: {line}")
            raise e

    def _init_dict(self):
        for sample_id in self.sample_id_list:
            self.match_dict[sample_id] = []

    def _approx_match(self, seq) -> str:
        return super()._approx_match(seq, self.fm, self.allow_mismatch)

    def _chunk_match_index(self, chunk, q):
        match_dict = {}
        undetermined = []
        for seq_id, seq in chunk:
            test_seq = seq[: self.barcode_len + len(self.pr)]
            all_poss_pr = [test_seq[i:i+len(self.pr)] for i in range(len(test_seq)-len(self.pr)+1)]
            fm_result = process.extractOne(self.pr, all_poss_pr, scorer=Levenshtein.distance, score_cutoff=int(len(self.pr)*0.15))
            if fm_result:
                seq = seq[: fm_result[2]]
            else:
                undetermined.append(seq_id)
                continue

            match_sample_id = self._approx_match(seq)

            if match_sample_id:
                if match_sample_id not in match_dict:
                    match_dict[match_sample_id] = []
                match_dict[match_sample_id].append(seq_id)
            else:
                undetermined.append(seq_id)

        q.put((match_dict, undetermined))

    def _match_index(self, fq_gz):
        self._get_chunks(fq_gz)

        queue = mp.Queue()
        processes = []

        for chunk in self.chunks:
            p = mp.Process(target=self._chunk_match_index, args=(chunk, queue,))
            p.start()
            processes.append(p)

        for _ in tqdm(range(len(self.chunks)), desc="Matching Index..."):
            match_dict, undetermined = queue.get()

            for sample_id, seq_ids in match_dict.items():
                self.match_dict[sample_id].extend(seq_ids)
            self.undetermined.extend(undetermined)

            del match_dict, undetermined

        for p in processes:
            p.join()

    def _rm_empty(self):
        self.match_dict = super()._rm_empty(self.match_dict)

    def _get_records(self, fq_gz):
        self.records = super()._get_records(fq_gz, fq_gz.replace(".gz", ""))

    def _write_match_fq(self):
        for sample_id, seq_id_list in self.match_dict.items():
            save_path = f"{self.save_dir}/{sample_id}_{self.suffix}.fastq.gz"
            super()._write_fq(self.records, seq_id_list, save_path)

    def _report(self):
        match_num = 0
        undetermined_num = len(self.undetermined)

        logger.info("Sample_ID\tCounts")
        for sample_id, seq_id_list in self.match_dict.items():
            logger.info(f"{sample_id}\t{len(seq_id_list)}")
            match_num += len(seq_id_list)

        if match_num+undetermined_num != 0:
            logger.info(f"Undetermined\t{undetermined_num}")
            logger.info(f"Determined rate: {match_num / (match_num + undetermined_num) * 100:.2f}%")
        else:
            logger.info("No reads matched.")

    def run(self, fq_gz, dry):
        self._match_index(fq_gz)
        self._rm_empty()
        self._get_records(fq_gz)
        if not dry:
            self._write_match_fq()
            super()._write_fq(self.records, self.undetermined, f"{self.save_dir}/Undetermined_{self.suffix}.fastq.gz")
        self._report()


class PairedEndDemux(Demux):
    def __init__(self,
            index_path,
            for_raw_path: str,
            rev_raw_path: str,
            save_dir: str,
            pr_5: str = "GTCGGTAAAACTCGTGCCAGC",
            pr_3: str = "CATAGTGGGGTATCTAATCCCAGTTTG",
            allow_mismatch: int = 0,
            threads: int = None,
            dry: bool = False
        ):
        """
        Demultiplex pair-end fastq.gz file.

        :param index_path: path to index file.
        :param for_raw_path: path to raw _R1.fastq.gz file.
        :param rev_raw_path: path to raw _R2.fastq.gz file.
        :param save_dir: The directory to save demultiplexed fastq.gz file.
        :param pr_5: 5' primer after forward index. Default is "GTCGGTAAAACTCGTGCCAGC"(MiFish-UF).
        :param pr_3: 3' primer after reverse index. Default is "CATAGTGGGGTATCTAATCCCAGTTTG"(MiFish-UR).
        :param allow_mismatch: allow distance for barcode matching. Default is 0.
        :param threads: Number of CPU cores to be used for processing. If not specify, use all. Default is None.
        :param dry: Dry run. If True, no files will be written; only a demultiplexed report will be output. Default is False.
        """
        super().__init__(save_dir, allow_mismatch, threads)
        self.match_dict = {
            "FF": {}, # seqs in forward file match forward index
            "RR": {}, # seqs in reverse file match reverse index
            "FR": {}, # seqs in forward file match reverse index
            "RF": {}, # seqs in reverse file match forward index
        }
        self.unmatch_dict = {"F": [], "R": []} # seqs did not matched  both in forward and reverse file
        self.undetermined_dict = {"F": [], "R": []} # seqs did not match any index in forward or reverse file
        self.match_num = {}
        self.unmatch_num = {}
        self._read_barcode(index_path)
        self._init_dict()
        self.for_fm = FuzzyMatch(self.for_barcode)
        self.rev_fm = FuzzyMatch(self.rev_barcode)
        self.pr = [pr_5, pr_3]
        self.run(for_raw_path, rev_raw_path, dry)

    def _read_barcode(self, barcode_path):
        try:
            self.sample_id_list = []
            self.for_barcode = {}
            self.rev_barcode = {}
            with open(barcode_path, 'r') as f:
                for line in tqdm(f, desc="Reading Index file"):
                    line = line.strip().split("\t")
                    self.sample_id_list.append(line[0])
                    self.for_barcode[line[0]] = line[1]
                    self.rev_barcode[line[0]] = line[2]
            self.barcode_len = len(line[1])
        except Exception as e:
            logger.error(f"Error reading Index file: {e}. Error line: {line}")
            raise e

    def _init_dict(self):
        for sample_id in self.sample_id_list:

            for m_dict in self.match_dict.values():
                m_dict[sample_id] = []

            self.match_num[sample_id] = 0
            self.unmatch_num[sample_id] = 0 

    def _approx_match(self, seq) -> str:
        for fm, direction in [(self.for_fm, "F"), (self.rev_fm, "R")]:
            match_sample_id = super()._approx_match(seq, fm, self.allow_mismatch)
            if match_sample_id:
                return match_sample_id, direction
        return None, None

    def _chunk_match_index(self, chunk, q, read_direct):
        match_dict = {
            "FF": {}, "RR": {}, 
            "FR": {}, "RF": {}, 
        }
        undetermined_dict = {"F": [], "R": []}

        for seq_id, seq in chunk:
            poss_pr_start = []
            for pr in self.pr:
                test_seq = seq[: self.barcode_len+len(pr)]
                all_poss_pr = [test_seq[i:i+len(pr)] for i in range(len(test_seq)-len(pr)+1)]
                fm_result = process.extractOne(pr, all_poss_pr, scorer=Levenshtein.distance, score_cutoff=int(len(pr)*0.15))
                if fm_result:
                    poss_pr_start.append(fm_result[2])
            if poss_pr_start != []:
                seq = seq[: min(poss_pr_start)]
            else:
                undetermined_dict[read_direct].append(seq_id)
                continue

            match_sample_id, index_direct = self._approx_match(seq)

            if match_sample_id:
                if match_sample_id not in match_dict[f"{read_direct}{index_direct}"]:
                    match_dict[f"{read_direct}{index_direct}"][match_sample_id] = []
                match_dict[f"{read_direct}{index_direct}"][match_sample_id].append(seq_id)
            else:
                undetermined_dict[read_direct].append(seq_id)

        q.put((match_dict, undetermined_dict))

    def _match_index(self, fq_gz, read_direct):
        self._get_chunks(fq_gz)

        queue = mp.Queue()

        processes = []
        for chunk in tqdm(self.chunks, desc="Assigning match tasks..."):
            p = mp.Process(target=self._chunk_match_index, args=(chunk, queue, read_direct,))
            p.start()
            processes.append(p)

        for _ in tqdm(range(len(self.chunks)), desc="Matching Index..."):
            match_dict, undetermined_dict = queue.get()

            for direct, match in match_dict.items():
                for sample_id, seq_ids in match.items():
                    self.match_dict[direct][sample_id].extend(seq_ids)
            for read_direct, seq_ids in undetermined_dict.items():
                self.undetermined_dict[read_direct].extend(seq_ids)

            del match_dict, undetermined_dict

        for p in processes:
            p.join()

    def _rm_empty(self):
        for key, match in self.match_dict.items():
            self.match_dict[key] = super()._rm_empty(match)

    def _find_matches(self, forward_matches, reverse_matches):
        common_sample_ids = set(forward_matches) & set(reverse_matches)
        self.common_match_dict = {}
        for sample_id in tqdm(common_sample_ids, desc="Searching Forward & Reverse common matches..."):
            common_seq_ids = list(set(forward_matches[sample_id]) & set(reverse_matches[sample_id]))
            self.common_match_dict[sample_id] = common_seq_ids
            self.match_num[sample_id] += len(common_seq_ids)

    def _find_unmatches(self, for_match_dict, rev_match_dict, for_unmatches, rev_unmatches):
        for sample_id, forward_seq_ids in tqdm(for_match_dict.items(), desc="Searching Forward unmatches..."):
            if sample_id not in rev_match_dict:
                for_unmatches.extend(forward_seq_ids)
                self.unmatch_num[sample_id] += len(forward_seq_ids)
            else:
                only_forward_seq_ids = list(set(forward_seq_ids) - set(rev_match_dict[sample_id]))
                for_unmatches.extend(only_forward_seq_ids)
                self.unmatch_num[sample_id] += len(only_forward_seq_ids)

        for sample_id, reverse_seq_ids in tqdm(rev_match_dict.items(), desc="Searching Reverse unmatches..."):
            if sample_id not in for_match_dict:
                rev_unmatches.extend(reverse_seq_ids)
                self.unmatch_num[sample_id] += len(reverse_seq_ids)
            else:
                only_reverse_seq_ids = list(set(reverse_seq_ids) - set(for_match_dict[sample_id]))
                rev_unmatches.extend(only_reverse_seq_ids)
                self.unmatch_num[sample_id] += len(only_reverse_seq_ids)

    def _get_records(self, for_fq_gz, rev_fq_gz):
        for_fq = for_fq_gz.replace(".gz", "")
        self.for_records = super()._get_records(for_fq_gz, for_fq)
        rev_fq = rev_fq_gz.replace(".gz", "")
        self.rev_records = super()._get_records(rev_fq_gz, rev_fq)

    def _write_match_fq(self, records, suffix):
        for sample_id, seq_id_list in self.common_match_dict.items():
            save_path = os.path.join(self.save_dir, f"{sample_id}_{suffix}.fastq.gz")
            super()._write_fq(records, seq_id_list, save_path)

    def _report(self):
        fh = logging.FileHandler(os.path.join(self.save_dir, "report.txt"))
        fh.setLevel(logging.INFO)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

        sample_stats = {}
        total_unmatch_num = 0
        undetermined_num = len(set(self.undetermined_dict["F"] + self.undetermined_dict["R"]))
        total_num = undetermined_num

        for sample_id in tqdm(self.sample_id_list, desc="Generating report..."):
            match_num = self.match_num[sample_id]
            unmatch_num = self.unmatch_num[sample_id]
            num = match_num + unmatch_num
            if num != 0:
                sample_stats[sample_id] = (num, match_num, match_num / num)
                total_num += num
                total_unmatch_num += unmatch_num

        total_match_num = total_num - total_unmatch_num - undetermined_num

        if total_num != 0:
            logger.info("Sample_ID\tReads\tMatched_Reads\tMatched_Rate")
            for sample_id, (total, matched, match_rate) in sample_stats.items():
                logger.info(f"{sample_id}\t{total}\t{matched}\t{match_rate*100:.2f}%")
            logger.info(f"Total_matched\tReads:{total_match_num}\tRate:{(total_match_num / total_num)*100:.2f}%")
            logger.info(f"Total_Unmatched\tReads:{total_unmatch_num}\tRate:{(total_unmatch_num / total_num)*100:.2f}%")
            logger.info(f"Total_Undetermined\tReads:{undetermined_num}\tRate:{(undetermined_num / total_num)*100:.2f}%")
        else:
            logger.info("No reads matched any index.")

    def run(self, for_fq_gz, rev_fq_gz, dry):
        self._match_index(for_fq_gz, "F")
        self._match_index(rev_fq_gz, "R")
        self._rm_empty()

        if not dry:
            self._get_records(for_fq_gz, rev_fq_gz)

        self._find_matches(self.match_dict["FF"], self.match_dict["RR"])
        self._find_unmatches(self.match_dict["FF"], self.match_dict["RR"], self.unmatch_dict["F"], self.unmatch_dict["R"])
        del self.match_dict["FF"], self.match_dict["RR"]
        if not dry:
            self._write_match_fq(self.for_records, "R1")
            self._write_match_fq(self.rev_records, "R2")

        self._find_matches(self.match_dict["RF"], self.match_dict["FR"])
        self._find_unmatches(self.match_dict["RF"], self.match_dict["FR"], self.unmatch_dict["R"], self.unmatch_dict["F"])
        del self.match_dict["RF"], self.match_dict["FR"]
        if not dry:
            self._write_match_fq(self.rev_records, "R1")
            self._write_match_fq(self.for_records, "R2")

        if not dry:
            super()._write_fq(self.for_records, set(self.unmatch_dict["F"]), f"{self.save_dir}/unmatched_R1.fastq.gz")
            super()._write_fq(self.rev_records, set(self.unmatch_dict["R"]), f"{self.save_dir}/unmatched_R2.fastq.gz")

            super()._write_fq(self.for_records, set(self.undetermined_dict["F"]), f"{self.save_dir}/undetermined_R1.fastq.gz")
            super()._write_fq(self.rev_records, set(self.undetermined_dict["R"]), f"{self.save_dir}/undetermined_R2.fastq.gz")

        self._report()