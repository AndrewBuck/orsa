#ifndef _ORSA_DEBUG_
#define _ORSA_DEBUG_

#include <orsaTBB/malloc.h>

#include <QMutex>

#include <cstdarg> // definition of std::va_list

#include <string>
#include <vector>

#include <time.h>
#include <unistd.h>

#if (_POSIX_TIMERS > 0) && defined(_POSIX_MONOTONIC_CLOCK)
#define __ADVANCED_TIMER__
#endif

namespace orsa {

    class Debug {
    public:
        enum DebugType {
            DT_DEBUG,
            DT_WARNING,
            DT_ERROR
        };
    protected:
        virtual const char * __head__(const DebugType);

    public:
        static Debug * instance();

    protected:
        Debug();
    public:
        virtual ~Debug();

    protected:
        class DebugMessage {
        public:
            DebugMessage() {
                headSet = false;
                bodySet = false;
            }

        public:
            void setHead(const std::string & s) {
                head = s;
                headSet = true;
            }
        public:
            const std::string & getHead() const {
                return head;
            }
        public:
            bool headIsSet() const {
                return headSet;
            }
        protected:
            std::string head;
            bool headSet;

        public:
            void setBody(const std::string & s) {
                body = s;
                bodySet = true;
            }
        public:
            const std::string & getBody() const {
                return body;
            }
        public:
            bool bodyIsSet() const {
                return bodySet;
            }
        protected:
            std::string body;
            bool bodySet;
        };

    protected:
        typedef std::vector<DebugMessage> QueueType;
        QueueType queue;

    public:
        void initTimer();
    protected:
#ifdef __ADVANCED_TIMER__
        bool     tp_init;
        timespec tp_zero;  // used by clock_gettime
#endif // __ADVANCED_TIMER__

    public:
        // virtual void head(const DebugType type, const char * file, const int line);
        virtual void head(const DebugType type, const char * file, const int line, const char * function);
    public:
        void trace(const char * fmt, ...); // do not make this one virtual
    protected:
        virtual void vtrace(const char * fmt, std::va_list list);

    public:
        void create();
    public:
        void flush();
    protected:
        virtual void printOut(const std::string &);

    private:
        // a mutex that provides access serialization to the queue
        QMutex mutex;

    protected:
        static Debug * _instance;
    };

    // don't use this macro directly
#define __ORSA_DEBUG_MACRO__(type, format, ...)                         \
    orsa::Debug::instance()->create();                                  \
        orsa::Debug::instance()->head(type,__FILE__,__LINE__,__func__);	\
        orsa::Debug::instance()->trace(format, ##__VA_ARGS__);          \
        orsa::Debug::instance()->flush();

    // user ready macros
#define ORSA_DEBUG(format, ...)   __ORSA_DEBUG_MACRO__(orsa::Debug::DT_DEBUG,   format, ##__VA_ARGS__)
#define ORSA_WARNING(format, ...) __ORSA_DEBUG_MACRO__(orsa::Debug::DT_WARNING, format, ##__VA_ARGS__)
#define ORSA_ERROR(format, ...)   __ORSA_DEBUG_MACRO__(orsa::Debug::DT_ERROR,   format, ##__VA_ARGS__)

} // namespace orsa

#endif // _ORSA_DEBUG_
