const MAX_WORKER: usize = 4;

use std::collections::VecDeque;
use std::sync::Arc;
use std::sync::Mutex;

/// A generic work queue for work elements which can be trivially copied.
/// Any producer of work can add elements and any worker can consume them.
/// WorkQueue derives Clone so that it can be distributed among threads.
#[derive(Clone)]
pub struct WorkQueue<T: Send + Copy> {
    inner: Arc<Mutex<VecDeque<T>>>,
}

impl<T: Send + Copy> WorkQueue<T> {
    fn new() -> Self {
        Self {
            inner: Arc::new(Mutex::new(VecDeque::new())),
        }
    }

    fn get_work(&self) -> Option<T> {
        let maybe_queue = self.inner.lock();
        if let Ok(mut queue) = maybe_queue {
            queue.pop_front()
        } else {
            // There's a problem with the mutex.
            panic!("WorkQueue::get_work() tried to lock a poisoned mutex");
        }
    }

    fn add_work(&self, work: T) -> usize {
        if let Ok(mut queue) = self.inner.lock() {
            queue.push_back(work);
            queue.len()
        } else {
            panic!("WorkQueue::add_work() tried to lock a poisoned mutex");
        }
    }
}

pub struct SyncFlagTx {
    inner: Arc<Mutex<bool>>,
}

impl SyncFlagTx {
    fn set(&mut self, state: bool) -> Result<(), ()> {
        if let Ok(mut v) = self.inner.lock() {
            *v = state;
            Ok(())
        } else {
            Err(())
        }
    }
}

#[derive(Clone)]
pub struct SyncFlagRx {
    inner: Arc<Mutex<bool>>,
}

impl SyncFlagRx {
    fn get(&self) -> Result<bool, ()> {
        if let Ok(v) = self.inner.lock() {
            // Deref the MutexGuard to get at the bool inside
            Ok(*v)
        } else {
            Err(())
        }
    }
}

pub fn new_syncflag(initial_state: bool) -> (SyncFlagTx, SyncFlagRx) {
    let state = Arc::new(Mutex::new(initial_state));
    let tx = SyncFlagTx {
        inner: state.clone(),
    };
    let rx = SyncFlagRx {
        inner: state.clone(),
    };

    return (tx, rx);
}

pub fn main<T, U, V, F, G>(todo: &[T], mut do_work: F, acc0: V, mut accumulate_results: G) -> V
where
    T: Copy + Send,
    U: Copy + Send,
    F: FnMut(T) -> U + Copy + Send + Sync,
    G: FnMut(U, V) -> V,
{
    let queue: WorkQueue<T> = WorkQueue::new();

    use std::sync::mpsc::channel;
    let (results_tx, results_rx) = channel();

    let (mut more_jobs_tx, more_jobs_rx) = new_syncflag(true);

    use std::thread;
    let mut threads = Vec::new();

    for thread_num in 0..MAX_WORKER {
        let thread_queue = queue.clone();
        let thread_results_tx = results_tx.clone();
        let thread_more_jobs_rx = more_jobs_rx.clone();
        let thread_do_work = do_work.clone();

        let handle = thread::spawn(move || {
            while thread_more_jobs_rx.get().unwrap() {
                if let Some(work) = thread_queue.get_work() {
                    let result = thread_do_work(work);
                    match thread_results_tx.send((work, result)) {
                        Ok(_) => (),
                        Err(_) => {
                            break;
                        }
                    }
                }
                std::thread::yield_now();
            }
        });
        threads.push(handle);
    }

    let mut jobs_total = todo.len();
    for work in todo.iter() {
        queue.add_work(*work);
    }

    let mut res = acc0;
    while jobs_total > 0 {
        match results_rx.recv() {
            Ok(work_res) => {
                jobs_total -= 1;
                res = accumulate_results(work_res.1, res);
            }
            Err(_) => {
                panic!("All workers died unexpectedly.");
            }
        }
    }

    more_jobs_tx.set(false).unwrap();

    for handle in threads {
        handle.join().unwrap();
    }

    res
}
