import subprocess
import time

def transform_int2knumber(knint):
    knstr = str(knint)
    while len(knstr) < 5:
        knstr = "0" + knstr
    return "K" + knstr

knumbers = [transform_int2knumber(k) for k in range(1, 27110)]

start = time.time()

for knumber in knumbers[2000:27110]:
    url = "https://rest.kegg.jp/get/{}".format(knumber)
    outdir = "Data/KO/{}.txt".format(knumber)
    try:
        txt_out = subprocess.run(
                ['wget', url, '-O', outdir],
                text=True,
                stdout=subprocess.PIPE,
                check=True)
    except:
        #log_notKnumber = subprocess.run("echo {} >> not_knumber.log".format(knumber))
        f = open("not_knumber.log", "a")
        f.write(knumber + "\n")
        f.close() 
        
    wait = subprocess.run(
             "sleep .1s",
             shell = True)
    
end = time.time()
print(end - start)

def remove_empty_files(logfile, data_dir):
    with open(logfile, "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            path_delete = data_dir + line + ".txt"
            try:
                subprocess.run(["rm", path_delete])
            except:
                continue
remove_empty_files("not_knumber.log", "Data/KO/")