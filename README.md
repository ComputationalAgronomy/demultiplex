# Demux
A python script to demultiplex illumina reads.

### Clone The Repository
```
git clone https://github.com/ComputationalAgronomy/demultiplex.git
```

## Installation
1. Required Python Package Installation
```sh
python -m pip install -r requirements.txt
```
2. Local package installs
```sh
python -m pip install -e .
```

## Usage

### Single-end
```python
from demux import SingleEndDemux

demux = SingleEndDemux(
    index_path = "/path/to/index.txt",
    raw_path = "/path/to/rawdata.fastq.gz",
    save_dir = "/path/to/save/dir/",
    allow_mismatch = 0,
    read_direction = "forward"
)
```
`index_path`: path to index .TXT file, format: `SAMPLE_ID<tab>INDEX`

`raw_path`: path to raw .FASTQ.GZ file.

`save_dir`: The directory to save demultiplexed .FASTQ.GZ file.

`allow_mismatch`: allow mismatch for barcode matching. Default is `0`.

`read_direction`: `forward` or `reverse`. Default is `"forward"`

### Paired-end
```python
from demux import PairedEndDemux

demux = PairedEndDemux(
    index_path = "/path/to/index.txt",
    for_raw_path = "/path/to/rawdata_R1.fastq.gz",
    rev_raw_path = "/path/to/rawdata_R2.fastq.gz",
    save_dir = "/path/to/save/dir/",
    allow_mismatch = 0
)
```
`index_path`: path to index .TXT file, format: `SAMPLE_ID<tab>FORWARD_INDEX<tab>REVERSE_INDEX`

`for_raw_path`: path to raw _R1.FASTQ.GZ file.

`rev_raw_path`: path to raw _R2.FASTQ.GZ file.

### Example
You can use the `example.ipynb` notebook and the materials in the `example` folder to practice.