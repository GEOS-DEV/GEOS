/*
 * exceptiontrace.cpp
 *
 *  Created on: Nov 22, 2014
 *      Author: rrsettgast
 */




#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#ifndef __USE_GNU
#define __USE_GNU
#endif

#include <execinfo.h>
#include <signal.h>
#include <string.h>

#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <sys/ucontext.h>

void my_terminate(void);

namespace {
    // invoke set_terminate as part of global constant initialization
    static const bool SET_TERMINATE = std::set_terminate(my_terminate);
}

/*
// This structure mirrors the one found in /usr/include/asm/ucontext.h
typedef struct _sig_ucontext {
   unsigned long     uc_flags;
   struct ucontext   *uc_link;
   stack_t           uc_stack;
   struct sigcontext uc_mcontext;
   sigset_t          uc_sigmask;
} sig_ucontext_t;
*/

void crit_err_hdlr(int sig_num, siginfo_t * info, void * ucontext_inp) {

  //sig_ucontext_t * uc = (sig_ucontext_t *)ucontext;
  __darwin_ucontext64 * uc = (__darwin_ucontext64 *)ucontext_inp;

    // Get the address at the time the signal was raised from the EIP (x86)
    void * caller_address = (void *) uc->uc_mcontext.eip;

    std::cerr << "signal " << sig_num
              << " (" << strsignal(sig_num) << "), address is "
              << info->si_addr << " from "
              << caller_address << std::endl;

    void * array[50];
    int size = backtrace(array, 50);

    std::cerr << __FUNCTION__ << " backtrace returned "
              << size << " frames\n\n";

    // overwrite sigaction with caller's address
    array[1] = caller_address;

    char ** messages = backtrace_symbols(array, size);

    // skip first stack frame (points here)
    for (int i = 1; i < size && messages != NULL; ++i) {
        std::cerr << "[bt]: (" << i << ") " << messages[i] << std::endl;
    }
    std::cerr << std::endl;

    free(messages);

    exit(EXIT_FAILURE);
}

void my_terminate() {
    static bool tried_throw = false;

    try {
        // try once to re-throw currently active exception
        if (!tried_throw++) throw;
    }
    catch (const std::exception &e) {
        std::cerr << __FUNCTION__ << " caught unhandled exception. what(): "
                  << e.what() << std::endl;
    }
    catch (...) {
        std::cerr << __FUNCTION__ << " caught unknown/unhandled exception."
                  << std::endl;
    }

    void * array[50];
    int size = backtrace(array, 50);

    std::cerr << __FUNCTION__ << " backtrace returned "
              << size << " frames\n\n";

    char ** messages = backtrace_symbols(array, size);

    for (int i = 0; i < size && messages != NULL; ++i) {
        std::cerr << "[bt]: (" << i << ") " << messages[i] << std::endl;
    }
    std::cerr << std::endl;

    free(messages);

    abort();
}

int throw_exception() {
    // throw an unhandled runtime error
    throw std::runtime_error("RUNTIME ERROR!");
    return 0;
}

int foo2() {
    throw_exception();
    return 0;
}

int foo1() {
    foo2();
    return 0;
}

int main(int argc, char ** argv) {
    struct sigaction sigact;

    sigact.sa_sigaction = crit_err_hdlr;
    sigact.sa_flags = SA_RESTART | SA_SIGINFO;

    if (sigaction(SIGABRT, &sigact, (struct sigaction *)NULL) != 0) {
        std::cerr << "error setting handler for signal " << SIGABRT
                  << " (" << strsignal(SIGABRT) << ")\n";
        exit(EXIT_FAILURE);
    }

    foo1();

    exit(EXIT_SUCCESS);
}
