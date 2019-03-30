import os
import subprocess


#######
# SCRATCH

# all_files = list(filter(lambda x: os.path.isfile(x), os.listdir()))


# for f in all_files:
#     subprocess.call(["rclone","copy", f , "onedrive-abhi:ena-genomes", "-vv"])

#rclone copy ERR036201_1.fastq.gz onedrive-abhi:ena-genomes -vv


#######



#rclone copy ERR036201_1.fastq.gz onedrive-abhi:ena-genomes -vv

def has_fastq_in_name(string):
    if (string.find("fastq") == -1):
        #print("NO")
        return 0
    else:
        #print("YES")
        return 1



all_files = list(filter(lambda x: os.path.isfile(x), os.listdir()))
all_fastq_files = list(filter(lambda x:has_fastq_in_name(x), all_files))




all_files = os.listdir()
all_fastq_files = list(filter(lambda x:has_fastq_in_name(x), all_files))
all_aria2_files = list(filter(lambda x:has_aria2_in_name(x), all_files))
all_relevant_files = list(set(all_fastq_files) - set(all_aria2_files))



all_1gz_files = []
all_gz_files = []

for x in all_relevant_files:
    x1 = x.split(".")[-2]
    if x1 == "1":
        all_1gz_files.append(x)
    else:
        all_gz_files.append(x)


def move_file_gz_to_gz_files_folder(f):
    shutil.move(f,"./gz_files/")





for f in all_gz_files:
    #move_file_gz_to_gz_files_folder(f)
    print(f)
