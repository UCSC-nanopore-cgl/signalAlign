#!/usr/bin/env python
from __future__ import print_function
from multiprocessing import Manager, Process, current_process
from pathos.multiprocessing import ProcessingPool as Pool
from pathos.pools import _ProcessPool
from timeit import default_timer as timer

import random
import time

TOTAL_KEY = "total"
FAILURE_KEY = "failure"
MEM_USAGE_KEY = "mem_usage"

# here is an example of a class which will "perform work".  this could also just be a function that is invoked
# instead of an object that's created and then had a function run
class ExampleWork(object):
    def __init__(self, argument, iterable, failure_rate=None):
        self.argument = argument
        self.iterable = iterable
        self.failure_rate = failure_rate

    def work(self):
        if self.failure_rate is not None and random.uniform(0,1) < self.failure_rate:
            raise Exception("Simulated Failure")
        print("Success: {argument:%s, iterable:%s, failure_rate:%s}" % (self.argument, self.iterable, self.failure_rate))


# here is an example of a service.  it needs to iterate over the work queue (until the 'STOP' element is reached)
# the exception handling helps when issues come up.  it returns information to the caller via the done_queue.
# all services must have arguments work_queue and done_queue respectively; all custom args must be by keyword
def example_service(work_queue, done_queue, service_name="example_service"):
    # prep
    total_handled = 0
    failure_count = 0

    #catch overall exceptions
    try:
        for f in iter(work_queue.get, 'STOP'):
            # catch exceptions on each element
            try:
                # logging
                print("[{}] '{}' processing {}".format(service_name, current_process().name, f))

                # put your work here
                work = ExampleWork(**f)
                work.work()

            except Exception as e:
                # get error and log it
                message = e.message if e.message is not None and len(e.message) != 0 else "{}:{}".format(type(e), str(e))
                error = "{} '{}' failed with: {}".format(service_name, current_process().name, message)
                print("[{}] ".format(service_name) + error)
                done_queue.put(error)
                failure_count += 1

            # increment total handling
            total_handled += 1

    except Exception as e:
        # get error and log it
        message = e.message if e.message is not None and len(e.message) != 0 else "{}:{}".format(type(e), str(e))
        error = "{} '{}' critically failed with: {}".format(service_name, current_process().name, message)
        print("[{}] ".format(service_name) + error)
        done_queue.put(error)

    finally:
        # logging and final reporting
        print("[%s] '%s' completed %d calls with %d failures"
              % (service_name, current_process().name, total_handled, failure_count))
        done_queue.put("{}:{}".format(TOTAL_KEY, total_handled))
        done_queue.put("{}:{}".format(FAILURE_KEY, failure_count))


# the run service function takes in a service (function) with an iterable list, arguments to be passed to
# the "work" function done on each element, and an argument to append to these 'work' arguments to describe the element.
# it accepts an optional service_argument parameter which is passed directly to the service
def run_service(service, iterable, iterable_arguments, iterable_argument_name, worker_count,
                service_arguments={}, log_function=print):

    start = time.time()
    args = list(iterable_arguments.keys())
    args.append(iterable_argument_name)
    if log_function is not None:
        log_function("[run_service] running service {} with {} workers".format(service, worker_count))

    # setup workers for multiprocessing
    work_queue = Manager().Queue()
    done_queue = Manager().Queue()

    # add everything to work queue
    jobs = []
    for x in iterable:
        args = dict({iterable_argument_name: x},
                    **iterable_arguments)
        work_queue.put(args)

    # start workers

    for w in range(worker_count):
        p = Process(target=service, args=(work_queue, done_queue), kwargs=service_arguments)
        p.start()
        jobs.append(p)
        work_queue.put('STOP')

    # wait for threads to finish, then stop the done queue
    for p in jobs:
        p.join()
    done_queue.put('STOP')

    # if example service model is used, metrics can be gathered in this way
    messages = []
    total = 0
    failure = 0
    for f in iter(done_queue.get, 'STOP'):
        if f.startswith(TOTAL_KEY): total += int(f.split(":")[1])
        elif f.startswith(FAILURE_KEY): failure += int(f.split(":")[1])
        else: messages.append(f)

    # if we should be logging and if there is material to be logged
    if log_function is not None and (total + failure + len(messages)) > 0:
        log_function("[run_service] Summary {}:\n[run_service]\tTime: {}s\n[run_service]\tTotal: {}\n[run_service]\tFailure: {}"
                     .format(service, int(time.time() - start), total, failure))
        log_function("[run_service]\tMessages:\n[run_service]\t\t{}".format("\n[run_service]\t\t".join(messages)))

    # return relevant info
    return total, failure, messages, None


def run_service2(service, iterable, iterable_arguments, iterable_argument_names, worker_count,
                 service_arguments={}, log_function=print):

    start = time.time()
    args = list(iterable_arguments.keys())
    args.extend(iterable_argument_names)
    if log_function is not None:
        log_function("[run_service] running service {} with {} workers".format(service, worker_count))

    # setup workers for multiprocessing
    work_queue = Manager().Queue()
    done_queue = Manager().Queue()

    # add everything to work queue
    jobs = []
    for x in iterable:
        if type(x) is not tuple:
            x = [x]
        args = dict(dict(zip(iterable_argument_names, x)),
                    **iterable_arguments)
        work_queue.put(args)


    # start workers
    for w in range(worker_count):
        p = Process(target=service, args=(work_queue, done_queue), kwargs=service_arguments)
        p.start()
        jobs.append(p)
        work_queue.put('STOP')

    # wait for threads to finish, then stop the done queue
    for p in jobs:
        p.join()
    done_queue.put('STOP')

    # if example service model is used, metrics can be gathered in this way
    messages = []
    output = []
    total = 0
    failure = 0
    no_stop = True
    while no_stop:
        f = done_queue.get()
        if isinstance(f, str):
            if f == "STOP":
                no_stop = False
            elif f.startswith(TOTAL_KEY):
                total += int(f.split(":")[1])
            elif f.startswith(FAILURE_KEY):
                failure += int(f.split(":")[1])
            else:
                messages.append(f)
        else:
            output.append(f)
    # if we should be logging and if there is material to be logged
    if log_function is not None and (total + failure + len(messages)) > 0:
        log_function("[run_service] Summary {}:\n[run_service]\tTime: {}s\n[run_service]\tTotal: {}\n[run_service]\tFailure: {}"
                     .format(service, int(time.time() - start), total, failure))
        log_function("[run_service]\tMessages:\n[run_service]\t\t{}".format("\n[run_service]\t\t".join(messages)))

    # return relevant info
    return total, failure, messages, output


def run_service3(service, iterable, iterable_arguments, iterable_argument_names, worker_count,
                 log_function=print):

    start = timer()
    args = list(iterable_arguments.keys())
    args.extend(iterable_argument_names)
    if log_function is not None:
        log_function("[run_service] running service {} with {} workers".format(service, worker_count))

    # add everything to work queue
    all_args = []
    for x in iterable:
        if type(x) is not tuple:
            x = [x]
        args = dict(dict(zip(iterable_argument_names, x)),
                    **iterable_arguments)
        all_args.append(args)

    pool = Pool(worker_count)
    results = pool.amap(service, all_args)
    final_results = results.get()
    # if example service model is used, metrics can be gathered in this way
    messages = []
    total = len(final_results)
    failure = 0
    for error, mem_usage in final_results:
        if error is not False:
            failure += 1
            if type(error) is str:
                messages.append(error)

    # if we should be logging and if there is material to be logged
    if log_function is not None and (total + failure + len(messages)) > 0:
        log_function("[run_service] Summary {}:\n[run_service]\tTime: {}s\n[run_service]\tTotal: {}\n[run_service]\tFailure: {}"
                     .format(service, int(timer() - start), total, failure))
        log_function("[run_service]\tMessages:\n[run_service]\t\t{}".format("\n[run_service]\t\t".join(messages)))

    # return relevant info
    return total, failure, messages



# example of how to run service
def main():
    # no simulated error rate
    print("NO ERROR")
    run_service(example_service, range(12), {"argument": "12 count, no error, 3 workers"}, "iterable", 3)


    # half simulated error rate
    print("\nHALF ERROR")
    run_service(example_service, range(12), {"argument": "12 count, .5 error, 3 workers", "failure_rate": .5},
                "iterable", 3, {"service_name":"custom_name"})


    # half simulated error rate
    print("\nFULL ERROR")
    run_service(example_service, range(12), {"argument": "12 count, 1 error, 12 workers", "failure_rate": 1},
                "iterable", 12, {"service_name":"custom_name"})


if __name__ == "__main__":
    main()