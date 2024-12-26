#pragma once

#include <chrono>

class Timer {
   private:
    chrono::high_resolution_clock::time_point __StartTime, __LastTime,
        __EndTime;
    const char *__processName;

   public:
    Timer() {
        __StartTime = chrono::high_resolution_clock::now();
        __LastTime = __StartTime, __EndTime = __StartTime;
        __processName = "Unnamed";
    }

    explicit Timer(const char *processName) {
        __StartTime = chrono::high_resolution_clock::now();
        __LastTime = __StartTime, __EndTime = __StartTime;
        __processName = processName;
    }

    /// Refresh time
    void refresh_time() {
        __StartTime = chrono::high_resolution_clock::now();
        __LastTime = __StartTime, __EndTime = __StartTime;
    }

    /// Record current time
    void record_current_time() {
        __LastTime = __EndTime;
        __EndTime = chrono::high_resolution_clock::now();
    };

    /// Get the operation time from the last recorded moment
    double get_operation_time() {
        this->record_current_time();
        chrono::duration<double> elapsed = __EndTime - __LastTime;
        return elapsed.count();
    }

    /// Print the operation time from the last recorded moment
    void log_operation_time() {
        const double duration = this->get_operation_time();
        cout << "->Time used (sec): " << duration << '\n';
    }

    /// Print the operation time from the last recorded moment with a given name
    void log_operation_time(const char *operationName) {
        const double duration = this->get_operation_time();
        cout << "->Time used (sec) for operation [" << operationName
             << "]: " << duration << '\n';
    }

    /// Get the total time from the beginning
    double get_total_time() {
        this->record_current_time();
        chrono::duration<double> elapsed = __EndTime - __StartTime;
        return elapsed.count();
    }

    /// Print the total time from the beginning
    void log_total_time() {
        const double duration = this->get_total_time();
        cout << "--->Time used (sec) for process [" << __processName
             << "]: " << duration << '\n';
    }

    /// Print the total time from the beginning to the last recorded moment
    void log_sub_total_time() const {
        chrono::duration<double> elapsed = __EndTime - __StartTime;
        cout << "--->Time used (sec) for process [" << __processName
             << "]: " << elapsed.count() << '\n';
    }
};
