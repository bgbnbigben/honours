#ifndef __TSQUEUE_H__
#define __TSQUEUE_H__

#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <memory>
#include <vector>
 
template <class T>
class tsqueue {

private:
    std::queue<T> q_;
    mutable std::mutex mutex_;
    std::condition_variable emptyLock_;

public:
    void push(T);
    std::shared_ptr<T> pop();
    std::shared_ptr<T> tryPop();
    const bool empty() const;
};

template <class T>
inline void tsqueue<T>::push(T pValue) {
    std::lock_guard<std::mutex> lock(this->mutex_);
    this->q_.push(pValue);
    this->emptyLock_.notify_one();
} 

template <class T>
inline std::shared_ptr<T> tsqueue<T>::pop() {
    std::unique_lock<std::mutex> lock(this->mutex_);
    this->emptyLock_.wait(lock,[this] { return !this->q_.empty();});
    // std::shared_ptr<T> ret(new T(this->q_.front()) performs 2 memory allocations!
    std::shared_ptr<T> ret=std::shared_ptr<T>(std::make_shared<T>(this->q_.front()));
    this->q_.pop();
    return ret;
}

template <class T>
inline std::shared_ptr<T> tsqueue<T>::tryPop() {
    if(this->q_.empty())
        return std::shared_ptr<T>(); // null
    std::shared_ptr<T> ret;
    {
        std::lock_guard<std::mutex> lock(this->mutex_);
        ret=std::shared_ptr<T>(std::make_shared<T>(this->q_.front()));
        this->q_.pop();
    }
    return ret;
}

template <class T>
inline const bool tsqueue<T>::empty() const {
    std::lock_guard<std::mutex> lock(this->mutex_);
    return this->q_.empty();
}

extern template class tsqueue<std::vector<double>>;
template void tsqueue<std::vector<double, std::allocator<double>>>::push(std::vector<double, std::allocator<double>> pValue);
template std::shared_ptr<std::vector<double, std::allocator<double>>> tsqueue<std::vector<double, std::allocator<double>>>::pop();
template std::shared_ptr<std::vector<double, std::allocator<double>>> tsqueue<std::vector<double, std::allocator<double>>>::tryPop();
template const bool tsqueue<std::vector<double, std::allocator<double>>>::empty() const;

# endif
