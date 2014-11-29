import os, sys
import multiprocessing as mp
import queue

def runner(q):
    while True:
        try:
            dirName = q.get(block = False)
            os.chdir("./" + dirName)
            print("running ", "./" + dirName)
            os.system("~/new_execute")
            os.chdir("../")
        except queue.Empty:
            print("Queue is empty")
            break


if len(sys.argv) < 2:
    print("Not enough arguments")
    exit(1)

print(sys.argv[1])
numThreads = int(sys.argv[1])
folders = mp.Queue()
threads = []

print("Threads number: ", numThreads)

for dirName in os.listdir("."):
    if not os.path.isdir(dirName): continue
    print(dirName)
    folders.put(dirName)

for i in range(numThreads):
    t = mp.Process(target = runner, args = (folders,))
    threads.append(t)
    t.start()

try:
    for t in threads:
        t.join()
    print("Done!")
except KeyboardInterrupt:
    print("\nTerminating...")
    for t in threads:
        t.terminate()
