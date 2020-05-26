#include <sstream>
#include "utils.h"

// Convert a length of time in seconds to
// a human readable string
std::string seconds_to_human(int secs)
{
    int mins = secs / 60;
    int hrs  = mins / 60;
    int days = hrs  / 24;
    hrs  -=   days * 24;
    mins -=  (days * 24 + hrs) * 60;
    secs -= ((days * 24 + hrs) * 60 + mins) * 60;
    std::stringstream ss;
    ss << days << "d " << hrs  << "h " 
       << mins << "m " << secs << "s";
    return ss.str();
}
