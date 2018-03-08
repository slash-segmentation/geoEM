#pragma once

#include <cstdio>
#include <CGAL/Timer.h>
#include <boost/config.hpp>
#include <boost/thread.hpp>


///////////////////////////////////////////////////////////////////////////////////////////////////
// A command-line progress indicator. It displays the progress and time so far along with the
// estimated time remaining. The updates to the console are done in a separate thread once every
// second. Using this should have essentially no impact on performance.
///////////////////////////////////////////////////////////////////////////////////////////////////

class Progress
{
    std::FILE* _out;
    bool _is_console;
    CGAL::Timer _timer;
    const char* _text;
    const size_t _total;
    const int _total_digits;
    boost::thread* _thread;
    size_t _cur;

    
    static bool _get_is_console(std::FILE* out);
    static int _count_digits(size_t x);
    void _thread_update() const;
    void _cleanup();

public:
    // Create a progress indicator with the given text and total number of time update() will be called.
    // The timer/updater is not started, you MUST call start().
    inline Progress(const char* text, size_t total, std::FILE* out = stderr) : _out(out), _is_console(_get_is_console(out)), _timer(), _text(text), _total(total), _total_digits(_count_digits(total)), _thread(nullptr), _cur(0) {}
    inline ~Progress() { this->_cleanup(); }

    // Checks if the timer is running
    inline bool is_running() const { return this->_timer.is_running(); }

    // Get the current amount of time that the the progress has been running
    inline double time() const { return this->_timer.time(); }

    // Start the progress indicator
    void start();

    // Update the value of the progress indicator (won't be actually updated until later)
    BOOST_FORCEINLINE void update() { ++this->_cur; }

    // Call to stop the progress updater
    void done();
};
