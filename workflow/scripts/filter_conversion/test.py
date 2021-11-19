
import subprocess
cmd1 = "python3 create_model.py -i FASTQ_dataset.csv -o initial -k 5 -n 3 -s 0.8"
cmd2 = "python3 make_predictions.py -id FASTQ_dataset.csv -o predictions.csv -m initial.joblib"

#subprocess.call(['bash', '-c', cmd1])
subprocess.call(['bash', '-c', cmd2])
