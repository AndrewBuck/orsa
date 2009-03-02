#ifndef _SSWG_
#define _SSWG_

#include <orsa/cache.h>
#include <orsa/double.h>

class OutputTemplateEntry {
 public:	
  OutputTemplateEntry(const std::string & name,
		      const mpz_class   & size,
		      const bool          doUpload,
		      const bool          isOptional) :	
    fileName(name),
    fileSize(size),
    upload(doUpload),
    optional(isOptional) { }
 public:
  orsa::Cache<std::string> fileName;
  orsa::Cache<mpz_class>   fileSize;
  orsa::Cache<bool>        upload;
  orsa::Cache<bool>        optional;
};

typedef std::vector<OutputTemplateEntry> OutputTemplateVector;

#endif // _SSWG_
