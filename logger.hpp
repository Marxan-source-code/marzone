#pragma once

#include <chrono>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include "common.hpp"


namespace marzone {
    using namespace std;

    class Logger {
    public:
    Logger() {}

    // Returns verbosity
    int GetVerbosity() { return verbosity; }

    // For perf debugging purposes.
    void ShowTimePassedMs(const chrono::steady_clock::time_point& startTime) {
        chrono::steady_clock::time_point currentTime = chrono::steady_clock::now();
        uint64_t mseconds_passed = chrono::duration_cast<std::chrono::microseconds>(currentTime - startTime).count();
        cout << mseconds_passed << "\n";
    }

    /* * * *  ShowTimePassed displays the time passed so far * * * * */
    void ShowTimePassed(const chrono::high_resolution_clock::time_point& startTime) {
        chrono::high_resolution_clock::time_point currentTime = chrono::high_resolution_clock::now();
        uint64_t seconds_passed = chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();

        string plural_h = seconds_passed/3600==1? " ": "s";
        string plural_m = seconds_passed/60==1? " ": "s";

        stringstream msg;
        msg << "Time passed so far is ";

        // Format a user friendly message
        if (seconds_passed >= 60*60)
        {
            msg << " " << seconds_passed/3600 << " hour" << plural_h << " "
                << (seconds_passed/60)%60 << " min" << plural_m << " "
                << seconds_passed%60 << "secs \n";
        }
        else {
            if (seconds_passed >= 60)
            {
                msg << " " << seconds_passed/60 << " min" << plural_m << "  and "
                    << seconds_passed%60 << " secs \n";
            }
            else
            {
                msg << seconds_passed << " secs \n";
            }
        }

        // print message
        cout << msg.str();
        if (savelog)
            savelog << msg.str();
    }

    void StartDebugTraceFile(string my_savelogname) {
        if (verbosity > 2)
            fdebugtrace.open(my_savelogname);
    }

    void SetLogFile(string my_savelogname) {
        savelog.open(my_savelogname);
        if (!savelog.is_open())
        {
            ShowErrorMessage("Error: Cannot save to log file " + my_savelogname + "\n");
        } /* open failed */
    }

    void AppendDebugTraceFile(string sMess)
    {
        if (fdebugtrace)
        {
            if (verbosity > 2)
            {
                fdebugtrace << sMess;
                fdebugtrace.flush(); // persist message.
            }
        }
    }

    void AppendLogFile(string sMess)
    {
        if (savelog)
            savelog << sMess;
    }

    // ShowGenProg displays a general progress message when verbosity > 1 
    void ShowGenProg(string sMess) {
        if (verbosity > 1)
        {
            cout << sMess;
            if (savelog)
                savelog << sMess;
        }    
    }

    /* ShowGenProgInfo displays a general progress with information
    message when verbosity > 5 */
    void ShowGenProgInfo(string sMess) {
        if (verbosity > 5)
        {
            cout << sMess;
            if (savelog)
                savelog << sMess;
        }   
    }

    /* ShowDetailedProgress shows detailed progress information
    message when verbosity > 3 */
    void ShowDetProg(string sMess)
    {
        if (verbosity > 3)
        {
            cout << sMess;
            if (savelog)
                savelog << sMess;
        }  
    }

    /* ShowProg displays fundamental progress information. Basic run summary */
    void ShowProg(string sMess) {
        cout << sMess;
        if (savelog)
            savelog << sMess;
    }

    /* ShowWarningMessage displays a warning message no matter what verbosity level */
    void ShowWarningMessage(string sMess)
    {
        cout << sMess;
        if (savelog)
            savelog << sMess;
    }

    /* ShowErrorMessage displays an error message. No matter what verbosity these are
    always displayed. The program is terminated.*/
    void ShowErrorMessage(string sMess)
    {
        ShowWarningMessage(sMess);
        CloseLogFile();
        CloseDebugFile();

        // exit with error.
        exit(EXIT_FAILURE);
    }

    void CloseLogFile()
    {
        savelog.close();
    }

    void CloseDebugFile() {
        if (fdebugtrace.is_open())
        {
            fdebugtrace.close();
        }
    }

    // verbosity
    int verbosity;

    ofstream savelog;
    ofstream fdebugtrace;

    };

} // namespace marzone