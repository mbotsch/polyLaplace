#pragma once

#include <thread>
#include <atomic>
#include <utility>
#include <functional>

template<int NumThreads>
void threadHelper(std::function<void(std::atomic<int>*)> fn)
{
    using namespace std;
    atomic<int> cnt(0);
    
    array<thread, NumThreads> threads;
    
    for(int i = 0; i < NumThreads; ++i)
        threads[i] = thread(fn, &cnt);
    
    for(int i = 0; i < NumThreads; ++i)
        threads[i].join();
}

template<int NumThreads>
void threadHelper(std::function<void(const int, const int)> fn, const int N)
{
    using namespace std;
    atomic<int> cnt(0);
    
    array<thread, NumThreads> threads;
    
    auto fun = [&fn, &cnt, N](const int j)
    {
        int i = 0;
        while((i = cnt++) < N)
        {
            fn(i, j);
        }
    };
    
    
    for(int i = 0; i < NumThreads; ++i)
        threads[i] = thread(fun, i);
    
    for(int i = 0; i < NumThreads; ++i)
        threads[i].join();
}

template<int NumThreads>
void threadHelper(std::function<void(const int)> fn, const int N)
{
    using namespace std;
    atomic<int> cnt(0);
    
    array<thread, NumThreads> threads;
    
    auto fun = [&fn, &cnt, N]()
    {
        int i = 0;
        while((i = cnt++) < N)
        {
            fn(i);
        }
    };
    
    
    for(int i = 0; i < NumThreads; ++i)
        threads[i] = thread(fun);
    
    for(int i = 0; i < NumThreads; ++i)
        threads[i].join();
}

