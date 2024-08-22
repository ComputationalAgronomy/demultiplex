import gzip
from Bio import SeqIO
from dict_trie import Trie
from functools import partial
import logging
import multiprocessing as mp
from tqdm import tqdm


FORMAT = "%(message)s"
logger = logging.getLogger()
logger.setLevel(logging.INFO)
formatter = logging.Formatter(FORMAT)

sh = logging.StreamHandler()
sh.setLevel(logging.INFO)
sh.setFormatter(formatter)
logger.addHandler(sh)

class SingleEndDemux():
    def __init__(self,
            index_path,
            raw_path,
            save_dir,
            allow_mismatch: int = 0,
            threads: int = None,
            read_direction: str = "forward",
            dry: bool = False
        ):
        """
        Demultiplex single-end fastq.gz file to different sample fastq.gz files.

        :param index_path: path to barcode file, format: <sample_id>\\t<barcode>.
        :param raw_path: path to raw .fastq.gz file.
        :param save_dir: The directory to save demultiplexed fastq.gz file.
        :param allow_mismatch: allow distance for barcode matching. Default is 0.
        :param threads: Number of CPU cores to be used for processing. If not specify, use all. Default is None.
        :param read_direction: forward or reverse. Default is "forward"
        :param dry: dry run
        """
        self.save_dir = save_dir
        self.suffix = "R1" if read_direction in ["forward", "Forward", "F", "R1"] else "R2"
        self.allow_dist = allow_mismatch
        self.match_dict = {}
        self.undetermined = []
        self.threads = threads if threads else mp.cpu_count()
        self._read_barcode(index_path)
        self._init_dict()
        self.trie = Trie(self.barcode)
        self.run(raw_path, dry)

    def _read_barcode(self, barcode_path):
        try:
            logger.info("Reading Index file...")
            self.sample_id_list = []
            self.barcode = []
            with open(barcode_path, 'r') as f:
                for line in f:
                    line = line.strip().split("\t")
                    self.sample_id_list.append(line[0])
                    self.barcode.append(line[1])
            self.barcode_len = len(self.barcode[0])
        except Exception as e:
            logger.error(f"Error reading Index file: {e}. Error line: {line}")
            raise e

    def _init_dict(self):
        for sample_id in self.sample_id_list:
            self.match_dict[sample_id] = []
    
    def _approximate_match(seq, trie, allow_dist, query_id_list) -> str:
        match_barcode = trie.best_levenshtein(seq, allow_dist)
        if match_barcode:
            match_index = list(trie).index(match_barcode)
            match_sample_id = query_id_list[match_index]
            return match_sample_id
        return None

    def _get_chunk_match_dict(chunk, trie, allow_dist, query_id_list):
        match_dict = {}
        undetermined = []
        for seq_id, seq in chunk:
            match_sample_id = SingleEndDemux._approximate_match(
                seq,
                trie,
                allow_dist,
                query_id_list,
            )

            if match_sample_id:
                if match_sample_id not in match_dict:
                    match_dict[match_sample_id] = []
                match_dict[match_sample_id].append(seq_id)
            else:
                undetermined.append(seq_id)

        return match_dict, undetermined

    def _get_match_dict(self, fastq_gz):
        records = []
        with gzip.open(fastq_gz, 'rt') as in_handle:
            for record in tqdm(SeqIO.parse(in_handle, "fastq"), desc="Reading raw data..."):
                seq_id = record.id
                seq = str(record.seq)[: self.barcode_len]
                records.append((seq_id, seq))

        chunk_size = len(records) // self.threads + 1
        chunks = [records[i:i+chunk_size] for i in range(0, len(records), chunk_size)]
        with mp.Pool() as pool:
            func = partial(
                SingleEndDemux._get_chunk_match_dict,
                trie = self.trie,
                allow_dist = self.allow_dist,
                query_id_list = self.sample_id_list,
            )
            results = list(tqdm(pool.map(func, chunks), total=len(chunks), desc="Matching Index..."))

        for match_dict, undetermined in results:
            for sample_id, seq_ids in match_dict.items():
                self.match_dict[sample_id].extend(seq_ids)

            self.undetermined.extend(undetermined)

    def _remove_empty(self):
        self.match_dict = {k: v for k, v in self.match_dict.items() if v!=[]}

    def _write_fastq(fastq_gz, seq_id_list, save_path):
        with gzip.open(fastq_gz, 'rt') as in_handle, gzip.open(save_path, 'at') as out_handle:
            for record in SeqIO.parse(in_handle, "fastq"):
                if record.id in seq_id_list:
                    SeqIO.write(record, out_handle, "fastq")

    def _write_chunk_match_fastq(chunk, raw, save_dir, suffix):
        for sample_id, seq_id_list in chunk:
            save_path = f"{save_dir}/{sample_id}_{suffix}.fastq.gz"
            SingleEndDemux._write_fastq(raw, seq_id_list, save_path)

    def _write_match_fastq(self, raw):
        id_set = [(sample_id, seq_id_list) for sample_id, seq_id_list in self.match_dict.items()]

        chunk_size = len(id_set) // self.threads + 1
        chunks = [id_set[i:i+chunk_size] for i in range(0, len(id_set), chunk_size)]
        with mp.Pool() as pool:
            func = partial(
                SingleEndDemux._write_chunk_match_fastq,
                raw = raw,
                save_dir = self.save_dir,
                suffix = self.suffix,
            )
            tqdm(pool.map(func, chunks), total=len(chunks), desc="Writing Matches FASTQ...")

    def _write_undetermined_fastq(self, raw):
        logger.info("Writing Undetermined FASTQ...")
        SingleEndDemux._write_fastq(raw, self.undetermined, f"{self.save_dir}/undetermined.fastq.gz")

    def _report_match_rate(self):
        match_num = 0
        undetermined_num = len(self.undetermined)

        logger.info("Sample_ID\tCounts")

        for sample_id, seq_id_list in self.match_dict.items():
            logger.info(f"{sample_id}\t{len(seq_id_list)}")

            match_num += len(seq_id_list)

        logger.info(f"Undetermined\t{undetermined_num}")
        logger.info(f"Determined rate: {match_num / (match_num + undetermined_num) * 100:.2f}%")

    def run(self, raw, dry):
        self._get_match_dict(raw)
        self._remove_empty()
        if not dry:
            self._write_match_fastq(raw)
            self._write_undetermined_fastq(raw)
        self._report_match_rate()

class PairedEndDemux():
    def __init__(self,
            index_path,
            for_raw_path: str,
            rev_raw_path: str,
            save_dir: str,
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
        :param allow_mismatch: allow distance for barcode matching
        :param threads: Number of CPU cores to be used for processing. If not specify, use all. Default is None.
        :param dry: dry run
        """
        self.save_dir = save_dir
        self.allow_dist = allow_mismatch
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
        self.threads = threads if threads else mp.cpu_count()
        self._read_barcode(index_path)
        self._init_dict()
        self.for_trie = Trie(self.for_barcode)
        self.rev_trie = Trie(self.rev_barcode)
        self.run(for_raw_path, rev_raw_path, dry)

    def _read_barcode(self, barcode_path):
        try:
            logger.info("Reading Index file...")
            self.sample_id_list = []
            self.for_barcode = []
            self.rev_barcode = []
            with open(barcode_path, 'r') as f:
                for line in f:
                    line = line.strip().split("\t")
                    self.sample_id_list.append(line[0])
                    self.for_barcode.append(line[1])
                    self.rev_barcode.append(line[2])
            self.barcode_len = len(self.for_barcode[0])
        except Exception as e:
            logger.error(f"Error reading Index file: {e}. Error line: {line}")
            raise e

    def _init_dict(self):
        for sample_id in self.sample_id_list:

            for m_dict in self.match_dict.values():
                m_dict[sample_id] = []

            self.match_num[sample_id] = 0
            self.unmatch_num[sample_id] = 0 

    def _approximate_match(seq, for_trie, rev_trie, allow_dist, query_id_list) -> str:
        for trie, direction in [(for_trie, "F"), (rev_trie, "R")]:
            match_barcode = trie.best_levenshtein(seq, allow_dist)
            if match_barcode:
                match_index = list(trie).index(match_barcode)
                match_sample_id = query_id_list[match_index]
                return match_sample_id, direction
        return None, None

    def _get_chunk_match_dict(
            chunk,
            for_trie,
            rev_trie,
            allow_dist,
            query_id_list,
            read_direct
        ):
        match_dict = {
            "FF": {}, 
            "RR": {}, 
            "FR": {}, 
            "RF": {}, 
        }
        undetermined_dict = {"F": [], "R": []}
        for seq_id, seq in chunk:
            match_sample_id, index_direct = PairedEndDemux._approximate_match(
                seq,
                for_trie,
                rev_trie,
                allow_dist,
                query_id_list,
            )

            if match_sample_id:
                if match_sample_id not in match_dict[f"{read_direct}{index_direct}"]:
                    match_dict[f"{read_direct}{index_direct}"][match_sample_id] = []
                match_dict[f"{read_direct}{index_direct}"][match_sample_id].append(seq_id)
            else:
                undetermined_dict[read_direct].append(seq_id)

        return match_dict, undetermined_dict

    def _get_match_dict(self, fastq_gz, read_direct):
        records = []
        with gzip.open(fastq_gz, 'rt') as in_handle:
            for record in tqdm(SeqIO.parse(in_handle, "fastq"), desc="Reading raw data..."):
                seq_id = record.id
                seq = str(record.seq)[: self.barcode_len]
                records.append((seq_id, seq))

        chunk_size = len(records) // self.threads + 1
        chunks = [records[i:i+chunk_size] for i in range(0, len(records), chunk_size)]

        with mp.Pool() as pool:
            func = partial(
                PairedEndDemux._get_chunk_match_dict,
                for_trie = self.for_trie,
                rev_trie = self.rev_trie,
                allow_dist = self.allow_dist,
                query_id_list = self.sample_id_list,
                read_direct = read_direct,
            )
            results = list(tqdm(pool.map(func, chunks), total=len(chunks), desc="Matching Index..."))

        for match_dict, undetermined_dict in results:
            for direct, match in match_dict.items():
                for sample_id, seq_ids in match.items():
                    self.match_dict[direct][sample_id].extend(seq_ids)

            for read_direct, seq_ids in undetermined_dict.items():
                self.undetermined_dict[read_direct].extend(seq_ids)

    def _remove_empty(self):
        for key, match in self.match_dict.items():
            self.match_dict[key] = {k: v for k, v in match.items() if v!=[]}

    def _find_common_matches(self, forward_matches, reverse_matches):
        common_sample_ids = set(forward_matches) & set(reverse_matches)
        self.common_match_dict = {}
        for sample_id in common_sample_ids:
            common_seq_ids = list(set(forward_matches[sample_id]) & set(reverse_matches[sample_id]))
            self.common_match_dict[sample_id] = common_seq_ids
            self.match_num[sample_id] += len(common_seq_ids)

    def _find_unmatches(self, for_match_dict, rev_match_dict, for_unmatches, rev_unmatches):
        for sample_id, forward_seq_ids in for_match_dict.items():
            if sample_id not in rev_match_dict:
                for_unmatches.extend(forward_seq_ids)
                self.unmatch_num[sample_id] += len(forward_seq_ids)
            else:
                only_forward_seq_ids = list(set(forward_seq_ids) - set(rev_match_dict[sample_id]))
                for_unmatches.extend(only_forward_seq_ids)
                self.unmatch_num[sample_id] += len(only_forward_seq_ids)

        for sample_id, reverse_seq_ids in rev_match_dict.items():
            if sample_id not in for_match_dict:
                rev_unmatches.extend(reverse_seq_ids)
                self.unmatch_num[sample_id] += len(reverse_seq_ids)
            else:
                only_reverse_seq_ids = list(set(reverse_seq_ids) - set(for_match_dict[sample_id]))
                rev_unmatches.extend(only_reverse_seq_ids)
                self.unmatch_num[sample_id] += len(only_reverse_seq_ids)

    def _write_fastq(fastq_gz, seq_id_list, save_path):
        with gzip.open(fastq_gz, 'rt') as in_handle, gzip.open(save_path, 'at') as out_handle:
            for record in SeqIO.parse(in_handle, "fastq"):
                if record.id in seq_id_list:
                    SeqIO.write(record, out_handle, "fastq")

    def _write_chunk_match_fastq(chunk, for_raw, rev_raw, save_dir):
        for sample_id, seq_id_list in chunk:
            save_path = f"{save_dir}/{sample_id}_R1.fastq.gz"
            PairedEndDemux._write_fastq(for_raw, seq_id_list, save_path)
            save_path = f"{save_dir}/{sample_id}_R2.fastq.gz"
            PairedEndDemux._write_fastq(rev_raw, seq_id_list, save_path)

    def _write_match_fastq(self, for_raw, rev_raw):
        id_set = [(sample_id, seq_id_list) for sample_id, seq_id_list in self.common_match_dict.items()]

        chunk_size = len(id_set) // self.threads + 1
        chunks = [id_set[i:i+chunk_size] for i in range(0, len(id_set), chunk_size)]

        with mp.Pool() as pool:
            func = partial(
                PairedEndDemux._write_chunk_match_fastq,
                for_raw = for_raw,
                rev_raw = rev_raw,
                save_dir = self.save_dir
            )
            tqdm(pool.map(func, chunks), total=len(chunks), desc="Writing Matches FASTQ...")

    def _write_unmatch_fastq(self, for_raw, rev_raw):
        logger.info("Writing Unmatched FASTQ...")
        PairedEndDemux._write_fastq(for_raw, self.unmatch_dict["F"], f"{self.save_dir}/unmatched_R1.fastq.gz")
        PairedEndDemux._write_fastq(rev_raw, self.unmatch_dict["R"], f"{self.save_dir}/unmatched_R2.fastq.gz")

    def _write_undetermined_fastq(self, for_raw, rev_raw):
        logger.info("Writing Undetermined FASTQ...")
        PairedEndDemux._write_fastq(for_raw, self.undetermined_dict["F"], f"{self.save_dir}/undetermined_R1.fastq.gz")
        PairedEndDemux._write_fastq(rev_raw, self.undetermined_dict["R"], f"{self.save_dir}/undetermined_R2.fastq.gz")

    def _calculate_sample_stats(self):
        self.sample_stats = {}
        self.total_unmatch_num = 0
        self.undetermined_num = len(set(self.undetermined_dict["F"] + self.undetermined_dict["R"]))
        total_num = self.undetermined_num

        for sample_id in self.sample_id_list:
            match_num = self.match_num[sample_id]
            unmatch_num = self.unmatch_num[sample_id]
            num = match_num + unmatch_num
            match_rate = match_num / num
            self.sample_stats[sample_id] = (num, match_num, match_rate)
            total_num += num
            self.total_unmatch_num += unmatch_num

        self.unmatched_rate = self.total_unmatch_num / total_num
        self.undetermined_rate = self.undetermined_num / total_num
    
    def _log_sample_stats(self):
        logger.info("Sample_ID\tReads\tMatched_Reads\tMatched_Rate")
        for sample_id, (total, matched, match_rate) in self.sample_stats.items():
            logger.info(f"{sample_id}\t{total}\t{matched}\t{match_rate:.2f}")
        logger.info(f"Total_Unmatched\tReads:{self.total_unmatch_num}\tRate:{self.unmatched_rate:.2f}")
        logger.info(f"Total_Undetermined\tReads:{self.undetermined_num}\tRate:{self.undetermined_rate:.2f}")

    def _report_match_rate(self):
        self._calculate_sample_stats()
        self._log_sample_stats()

    def run(self, for_raw, rev_raw, dry):
        self._get_match_dict(for_raw, "F")
        self._get_match_dict(rev_raw, "R")
        self._remove_empty()

        self._find_common_matches(self.match_dict["FF"], self.match_dict["RR"])
        self._find_unmatches(self.match_dict["FF"], self.match_dict["RR"], self.unmatch_dict["F"], self.unmatch_dict["R"])
        if not dry:
            self._write_match_fastq(for_raw, rev_raw)

        self._find_common_matches(self.match_dict["RF"], self.match_dict["FR"])
        self._find_unmatches(self.match_dict["RF"], self.match_dict["FR"], self.unmatch_dict["R"], self.unmatch_dict["F"])
        if not dry:
            self._write_match_fastq(rev_raw, for_raw)
            self._write_unmatch_fastq(for_raw, rev_raw)
            self._write_undetermined_fastq(for_raw, rev_raw)
        
        self._report_match_rate()