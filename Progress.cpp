#include "Progress.hpp"

#include <boost/chrono.hpp>

#ifdef _WIN32
#ifdef _MSC_VER
#define PRIuPTR "Iu"
#else
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#endif
#include <io.h>
#else
#define __STDC_FORMAT_MACROS
#include <unistd.h>
#include <inttypes.h>
#endif

int Progress::_count_digits(size_t x) {
	int answer = 0;
	for (; x != 0; x /= 10) { ++answer; }
	return answer;
}
bool Progress::_get_is_console(FILE* out) { return isatty(fileno(out)) != 0; }

void Progress::start()
{
	if (this->_is_console)
	{
		std::fprintf(this->_out, "%s %*u/%" PRIuPTR "\r", this->_text, this->_total_digits, 0, this->_total);
		if (!this->_thread) { this->_thread = new boost::thread(boost::bind(&Progress::_thread_update, this)); }
	}
	this->_timer.start();
	this->_cur = 1;
}
void Progress::_thread_update() const
{
	const boost::chrono::milliseconds one_sec(999);
	for (;;)
	{
		boost::this_thread::sleep_for(one_sec);
		double t = this->_timer.time(), rem = t*(this->_total-this->_cur)/this->_cur;
		unsigned int t_min = ((unsigned int)t)/60, rem_min = ((unsigned int)rem)/60;
		std::fprintf(this->_out, "%s %*" PRIuPTR "/%" PRIuPTR " in %u:%02u - %u:%04.1f remaining   \r", this->_text, this->_total_digits, this->_cur, this->_total, t_min, (unsigned int)t - t_min*60, rem_min, rem - rem_min*60);
	}
}
void Progress::done()
{
	this->_timer.stop();
	this->_cleanup();
	double t = this->_timer.time();
	unsigned int t_min = ((unsigned int)t)/60;
	std::fprintf(this->_out, this->_is_console ? ("%s %" PRIuPTR "/%" PRIuPTR " in %u:%04.1f                         \n") : ("%s %" PRIuPTR "/%" PRIuPTR " in %u:%04.1f\n"), this->_text, this->_total, this->_total, t_min, t - t_min*60);
}
void Progress::_cleanup()
{
	if (this->_thread)
	{
		this->_thread->interrupt();
		this->_thread->join();
		delete this->_thread;
		this->_thread = nullptr;
	}
}
