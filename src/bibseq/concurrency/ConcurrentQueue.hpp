#pragma once
/*
 * ConcurrentQueue.hpp
 *
 *  Created on: Jan 7, 2016
 */

/*
 Copyright 2014 Arjan van der Velde, vandervelde.ag [at] gmail
 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at
 http://www.apache.org/licenses/LICENSE-2.0
 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 */

#pragma once

#include <queue>
#include <mutex>
#include <condition_variable>

namespace njhseq {
namespace concurrent {

// concurrent multi-producer/multi-consumer queue
template<typename T>
class ConcurrentQueue {
private:
	std::queue<T> queue_; // queue
	std::mutex mtx_; // lock
	std::condition_variable conditionvar_; // condition variable
public:
	// push onto the queue

	void push(T const& data) {
		{ // GUARD
			std::unique_lock<std::mutex> lock(mtx_);
			queue_.push(data);
		}
		conditionvar_.notify_one(); // NOTIFY
	}

	// isEmpty?
	bool empty()  {
		std::unique_lock<std::mutex> lock(mtx_); // GUARD
		return queue_.empty();
	}

	// try non-blocking pop
	bool tryPop(T& value) {
		std::unique_lock<std::mutex> lock(mtx_); // GUARD
		if (queue_.empty()) {
			return false;
		}
		value = queue_.front();
		queue_.pop();
		return true;
	}

	// wait and pop
	void waitPop(T& value) {
		std::unique_lock<std::mutex> lock(mtx_); // GUARD
		while (queue_.empty()) {
			conditionvar_.wait(lock); // WAIT
		}
		value = queue_.front();
		queue_.pop();
	}
};

}  // namespace concurrent
}  // namespace njhseq
