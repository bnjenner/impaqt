#include <functional>
#include <cstdint>
#include <queue>

// Special thanks to https://github.com/embeddedartistry/embedded-resources/blob/master/examples/cpp/dispatch.cpp 
//	for the code inspiration :)

extern bool MAIN_THREAD;

//////////////////////////////////////
// Thread safe dispatch queue
class thread_queue {

	typedef std::function<void(void)> call;

	public:

		std::mutex mlock;
		std::vector<std::thread> threads;
		std::queue<call> call_queue;
		std::condition_variable cv;
		bool quit = false;

		thread_queue(size_t thread_cnt) : threads(thread_cnt) {	
			for (size_t i = 0; i < threads.size(); i++) {
				threads[i] = std::thread(&thread_queue::thread_handler, this);
			}
		}

		~thread_queue() {
			
			// creat lock and set quit to true
			std::unique_lock<std::mutex> lock(mlock);  
			quit = true;

			// unlock, notify all threads
			lock.unlock();
			cv.notify_all();

			for (auto &t: threads) {
				if (t.joinable()) {
					t.join();
				}
			}
			::MAIN_THREAD = true;
		}

		// dispatch (basically enqueue)
		void dispatch(const call &&job) {		
			std::unique_lock<std::mutex> lock(mlock); 	// create lock
			call_queue.push(std::move(job));		 	// enqueue
			lock.unlock(); 								// unlock
			cv.notify_one();							// notify conditional var
		}

		// thread handler
		void thread_handler(void) {

			std::unique_lock<std::mutex> lock(mlock); // create lock

			do {
				// wait for data (still need to think about this one)
				cv.wait(lock, [this]{
					return (call_queue.size() || quit);
				});

				// after wait, we own the lock
				if (!quit && call_queue.size()) {

					// create job
					auto job = std::move(call_queue.front());
					call_queue.pop();

					lock.unlock();	// unlock now that we're done messing with the queue
					job(); 			// call job
					lock.lock(); 	// lock 
				}
			} while (!quit);
		
		}

		bool finished(){ return (call_queue.empty()); }
		int size(){ return call_queue.size(); }
};
