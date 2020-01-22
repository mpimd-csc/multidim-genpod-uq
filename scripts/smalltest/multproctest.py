from multiprocessing import Process, Queue


def f(val, pos, queue):
    queue.put((val, pos))


if __name__ == '__main__':
    alist = [0, 0]
    pqueue = Queue()
    pone = Process(target=f, args=(1, 1, pqueue))
    ptwo = Process(target=f, args=(2, 0, pqueue))
    pone.start()
    ptwo.start()
    pone.join()
    ptwo.join()
    athng = pqueue.get()
    alist[athng[1]] = athng[0]
    athng = pqueue.get()
    alist[athng[1]] = athng[0]
    print(alist)
