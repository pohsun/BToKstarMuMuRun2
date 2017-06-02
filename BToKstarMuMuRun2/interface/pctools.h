#ifndef COMMON_PCTOOLS_H
#define COMMON_PCTOOLS_H

// basic IO
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <string>

namespace pctools {
    class PMessanger {
        /*
         * Modes (STDIO/FIO/SIO) are supported.
         * For STDIO mode,
         *      The output is dumped to /dev/tty.
         * For FIO   mode,
         *      The output is dumped to log file,
         * For SIO   mode,
         *      The output is FORCELY showed on screen.
         */

        public:
            PMessanger(std::string);
            ~PMessanger();

            // main working functions
                // controlled by user
            void out(std::string, int verbosityLevel, int mode);
            void out(std::string, int verbosityLevel);
            void out(std::string);
                // 'sprint' on screen, 'fprint' for FIO
            void fprint(std::string, int verbosityLevel);
            void fprint(std::string);
            void sprint(std::string, int verbosityLevel);
            void sprint(std::string);
            void switchRedirectStdio(bool);
            void switchRedirectStdio();

            // constant, enum
                // CIRTICAL, ERROR always on stdout and file
            enum modes{ STDIO, FIO, SIO, TEE, nMODES};
            enum verbosityLevels{ CRITICAL, ERROR, SILENT, WARNING, INFO, HINT, DEBUG, nLEVELS};
            const std::string verbosityTitle[nLEVELS] = {
                "CRITICAL: ",
                "ERROR   : ",
                "SILENT  : ", // Not really used.
                "WARNING : ",
                "INFO    : ", // Default
                "HINT    : ",
                "DEBUG   : "};

            // isA, hasA

            // setters and getters
            int setMode(int);
            int setVerbosity(int);

            int  getMode() const;
            int  getVerbosity() const;
            bool isRedirectStdio() const;
            FILE* getOfile() const;

        private:
            std::string ofName_;
            bool isRedirectStdio_ = false;
            int currentMode_    = STDIO;
            int verbosityLevel_ = INFO;
            FILE *ofPtr_ = 0;

    };
};

#endif
