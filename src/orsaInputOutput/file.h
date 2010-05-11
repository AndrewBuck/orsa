#ifndef _ORSA_INPUT_OUTPUT_FILE_
#define _ORSA_INPUT_OUTPUT_FILE_

#include <zlib.h>

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <cstdio>
#include <string>

#include <errno.h>
#include <string.h>

#include <orsa/cache.h>
#include <orsa/debug.h>

namespace orsaInputOutput {
  
    class File : public osg::Referenced  {
    
    protected:
        enum FileStatus {
            FS_CLOSE  = 0,
            FS_OPEN_R = 1,
            FS_OPEN_W = 2,
            FS_OPEN_A = 3
        };
    
    public:
        File() : Referenced(true) { }
      
    protected:
        virtual ~File() { }
    
    public:
        virtual bool openR() = 0;
        virtual bool openW() = 0;
        virtual bool openA() = 0;
    
    public:
        virtual int close() = 0;
    
    public:
        virtual void rewind() = 0;
    
    public:
        virtual bool gets(char * buffer,int length) = 0;
    
    public:
        virtual bool puts(const char * buffer) = 0;
    
    public:
        virtual int read(void * buffer, size_t size, size_t length) = 0;
    
    public:
        virtual int write(const void * buffer, size_t size, size_t length) = 0;
    
    public:
        virtual int seek(long offset, int whence) = 0;
    
    public:
        virtual int flush() = 0;
    
    public:
        const std::string & getFileName() const { 
            return _filename.getRef();
        }
    
    public:
        virtual void setFileName(const std::string & filename) {
            close();
            _filename.set(filename);
        }
    
    
    protected:
        orsa::Cache<std::string> _filename;
    
    protected:
        orsa::Cache<FileStatus>  _status;
    };
  
    /*** NOTE: keep the code into CompressedFile and PlainFile in sync! ***/
  
    class CompressedFile : public orsaInputOutput::File {
    
    public:
        CompressedFile() : File() {
            _file = Z_NULL;
            _status.set(FS_CLOSE);
        }
      
    protected:     
        ~CompressedFile() { 
            close();
        }
    
    private:
        void printError(const char * func) const {
            int errnum;
            const char * err = gzerror(_file,&errnum);
            ORSA_ERROR("func: %s   file [%s]: %s   [_file: %x]",
                       func,
                       _filename.getRef().c_str(),
                       err,
                       _file);
        }
    
    private:
        bool _open(const char * mode) {
            close();
            _file = gzopen(_filename.getRef().c_str(),mode);
            if (_file == Z_NULL) {
                printError(__func__);
                return false;
            } else {
                return true;
            } 
        }
    
    public:
        bool openR() {
            const bool success = _open("r");
            if (success) {
                _status.set(FS_OPEN_R);
            }
            return success;
        }
    public:
        bool openW() { 
            const bool success = _open("w");
            if (success) {
                _status.set(FS_OPEN_W);
            }
            return success;
        }
    public:
        bool openA() {
            const bool success = _open("a");
            if (success) {
                _status.set(FS_OPEN_A);
            }
            return success;
        }
    
    public:
        int close() {
            if (_status.getRef() != FS_CLOSE) {
                const int retVal = gzclose(_file);
                if (retVal == Z_OK) {
                    _status.set(FS_CLOSE);
                    _file = 0;
                    return 0;
                } else {
                    printError(__func__);
                    return -1;
                }
            } else {
                // file already closed
                return 0;
            }
        }
    
    public:
        void rewind() {
            if (gzrewind(_file) != 0) {
                printError(__func__);
            }
        }
    
    public:
        bool gets(char * buffer,int length) {
            if (gzgets(_file,buffer,length) == Z_NULL) {
                return false;
            } else {
                return true;
            }
        }	
    
    public:
        bool puts(const char * buffer) {
            return (gzputs(_file,buffer) > 0);
        }
    
    public:
        int read(void * buffer, size_t size, size_t length) {
            return gzread(_file,buffer,size*length);
        }
    
    public:
        int write(const void * buffer, size_t size, size_t length) {
            return gzwrite(_file,buffer,size*length);
        }
    
    public:
        int seek(long offset, int whence) {
            return gzseek(_file,offset,whence);
        }
    
    public:
        int flush() {
            return gzflush(_file,Z_FULL_FLUSH);
        }
    
    private:
        gzFile _file;
    };
  
    //! Plain file, i.e. not compressed.
    class PlainFile : public orsaInputOutput::File {
    
    public:
        PlainFile() : File() {
            _file = 0;
            _status.set(FS_CLOSE);
        }
      
    protected:     
        ~PlainFile() { 
            close();
        }
    
    private:
        bool _open(const char * mode) {
            close();
            _file = fopen(_filename.getRef().c_str(),mode);
            if (_file == 0) {
                ORSA_ERROR("cannot open file [%s]: %s",
                           _filename.getRef().c_str(),
                           strerror(errno));
                return false;
            } else {
                return true;
            } 
        }
    
    public:
        bool openR() {
            const bool success = _open("r");
            if (success) {
                _status.set(FS_OPEN_R);
            }
            return success;
        }
    public:
        bool openW() { 
            const bool success = _open("w");
            if (success) {
                _status.set(FS_OPEN_W);
            }
            return success;
        }
    public:
        bool openA() {
            const bool success = _open("a");
            if (success) {
                _status.set(FS_OPEN_A);
            }
            return success;
        }
    
    public:
        int close() {
            if (_status.getRef() != FS_CLOSE) {
                if (fclose(_file) == 0) {
                    _status.set(FS_CLOSE);
                    _file = 0;
                    return 0;
                } else {
                    ORSA_ERROR("cannot close file [%s]: %s",
                               _filename.getRef().c_str(),
                               strerror(errno));
                    return -1;
                }
            } else {
                // file already closed
                return 0;
            }
        }
    
    public:
        void rewind() {
            ::rewind(_file);
        }
    
    public:
        bool gets(char * buffer,int length) {
            return (fgets(buffer,length,_file) != 0);
        }	
    
    public:
        bool puts(const char * buffer) {
            return (fputs(buffer,_file) >= 0);
        }
    
    public:
        int read(void * buffer, size_t size, size_t length) {
            return fread(buffer,size,length,_file);
        }
    
    public:
        int write(const void * buffer, size_t size, size_t length) {
            return fwrite(buffer,size,length,_file);
        }
    
    public:
        int seek(long offset, int whence) {
            return fseek(_file,offset,whence);
        }
    
    public:
        int flush() {
            return fflush(_file);
        }
    
    private:
        FILE * _file;
    };
  
  
    //! T should be either a CompressedFile or a PlainFile
    //! D is the data type, for data storage, i.e. std::list<double> or std::vector<double>
    template <class T, class D> class InputFile : public osg::Referenced {
    
    public:
        InputFile() : osg::Referenced(true) {
            _lineLength.set(1024);
            _file = new T;
            dataInit();
        }
      
    protected:
        virtual ~InputFile() {
            dataDestroy();
        }
    
    public:
        const std::string & getFileName () const { 
            return _file->getFileName();
        }
    
    public:
        virtual void setFileName (const std::string & filename) {
            _file->setFileName(filename);
        }
    
    public:
        virtual bool read() {
            if (!_file->openR()) {
                return false;
            }
            const unsigned int length = _lineLength.getRef();
            char line[length], lineAbove[length];
            _file->rewind();
            while (_file->gets(line,length)) {
                line[strlen(line)-1] = '\0'; // remove trailing \n
                if (goodLine(line)) {
                    // uncomment while debugging
                    // ORSA_DEBUG("accepted line: [%s]",line);
                    // first, try to process a single line
                    // ORSA_DEBUG("calling processLine, line: [%s]",line);
                    if (!processLine(line)) {
                        // if unsuccessful, try processing two lines together
                        // this is necessary as sometimes files have mixed 1-line, 2-line data
                        if (goodLine(lineAbove)) {
                            // ORSA_DEBUG("calling processLines, line1: [%s]",lineAbove);
                            // ORSA_DEBUG("calling processLines, line2: [%s]",line);
                            if (!processLines(lineAbove,line)) {
                                // ORSA_DEBUG("processLines failed");
                            }
                        } else {
                            // ORSA_DEBUG("processLine failed");
                        }
                    }
                    //
                    strcpy(lineAbove,line);
                } else {
                    // uncomment while debugging
                    // ORSA_DEBUG("rejected line: [%s]",line);
                }
            }
            _file->close();
            return true;
        }
    
    protected:
        //! methods needed to allocate and free memory for _data, when needed;
        virtual void dataInit() { }
        virtual void dataDestroy() { }
    
    public:
        //! basic checks, i.e. is not a comment line
        virtual bool goodLine(const char * line) = 0;

    public:
        //! extract data from a good line
        virtual bool processLine(const char * line) = 0;

    public:
        //! extract data from two good lines,
        //! defaults to empty method returning false
        virtual bool processLines(const char * /* line_1 */ , const char * /* line_2 */ ) { return false; }
    
    public:
        orsa::Cache<unsigned int> _lineLength;
    private:
        osg::ref_ptr<T> _file;
    public:
        D _data;
    };
  
} // namespace orsaInputOutput

#endif // _ORSA_INPUT_OUTPUT_FILE_
