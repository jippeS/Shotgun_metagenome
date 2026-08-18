start_script=$SECONDS

# set tmp
mkdir -p /export/microlab/temp/temp_${SLURM_JOB_ID}
export TMPDIR="/export/microlab/temp/temp_${SLURM_JOB_ID}"
echo "Designated temporary directory for job $SLURM_JOB_ID: $TEMPDIR"

## activate conda environment

## On Linux 
echo "Activating Bedtools conda environment"
source /export/microlab/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-metagenome-2024.5

echo -e "Print QIIME2 version info to stdout"

# print q2 version info
qiime info

#----------------------------------------------------------------------------
## DEFINE VARIABLES

#DNA or RNA
type="DNA"
cpus=16
compile_results=no

# Define working directory for sample input and output
sample_dir="/export/microlab/users/PVEE/IDIN"

## DEFINE SAMPLES
# Validate type
if [ "$type" != "DNA" ] && [ "$type" != "RNA" ]; then
  echo "Invalid type specified. Please use 'DNA' or 'RNA'."
  exit 1
fi

# Validate cpus
if [ "$cpus" -lt 1 ] || [ "$cpus" -gt 16 ]; then
  echo "Invalid number of CPUs specified. Please use a value between 1 and 16."
  exit 1
fi

# Define the base directory for reads of either DNA or RNA sequences
if [ "$type" == "DNA" ]; then
  reads_path="/export/projects/Transcriptomic/DNA/Dataset1_23052024/EN00005451_hdd1/IDIN"
  read_length=150
  echo -e "Reads path set for type 'DNA metagenomic 150bp PE sequence data' to:\n$reads_path"
else
  reads_path="/export/projects/Transcriptomic/RNA/Dataset1_23052024/EN00005454_hdd1/IDIN"
  read_length=100
  echo -e "Reads path set for type 'RNA metatranscriptomic 100bp PE sequence data' to:\n$reads_path"
fi  


# set tempdir (first remove/empty any temp files present)
mkdir -p /export/microlab/temp/temp_${SLURM_JOB_ID}
export TMPDIR="/export/microlab/temp/temp_${SLURM_JOB_ID}"
echo "Designated temporary directory for job $SLURM_JOB_ID: $TEMPDIR"


# Read the DNA sample list into an array
sample_list=("IDIN02" "IDIN03" "IDIN06" "IDIN07" "IDIN11" "IDIN12" "IDIN18" "IDIN19" "IDIN28" "IDIN38" "IDIN39") #($(<"$sample_dir/${type}_samples.tsv")) #("IDIN02" "IDIN12" "IDIN19")

# create dir for merged data
mkdir -p "${sample_dir}/IDIN_${type}_merged"
merged_dir="${sample_dir}/IDIN_${type}_merged"

# create cache for merged data
if [ ! -d "${merged_dir}/cache" ]; then
      echo -e "Creating q2 cache for $sample..."

      qiime tools cache-create \
        --cache "${merged_dir}/cache"

else
    echo -e "q2 cache for merged data already exists"
fi 
