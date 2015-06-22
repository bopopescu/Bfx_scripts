import multiprocessing
from multiprocessing import Process, Queue
from time import sleep
from random import uniform

def doCalculation(par):
    t = uniform(0,2)
    sleep(t)
    return par * par  # just to simulate some calculation

def feed(queue, parlist):
    for par in parlist:
            queue.put(par)

def calc(queueIn, queueOut):
    while True:
        try:
            par = queueIn.get(block = False)
            print "dealing with ", par, "" 
            res = doCalculation(par)
            queueOut.put((par,res))
        except:
            break

def write(queue, fname):
    fhandle = open(fname, "w")
    while True:
        try:
            par, res = queue.get(block = False)
            print >>fhandle, par, res
        except:
            break
    fhandle.close()





if __name__ == "__main__":
    nthreads = multiprocessing.cpu_count()
    fname = "foo"
    workerQueue = Queue()
    writerQueue = Queue()
    parlist = [1,2,3,4,5,6,7,8,9,10]
    feedProc = Process(target = feed , args = (workerQueue, parlist))
    calcProc = [Process(target = calc , args = (workerQueue, writerQueue)) for i in range(nthreads)]
    writProc = Process(target = write, args = (writerQueue, fname))


    feedProc.start()
    feedProc.join ()
    for p in calcProc:
        p.start()
    for p in calcProc:
        p.join()
    writProc.start()
    writProc.join ()
