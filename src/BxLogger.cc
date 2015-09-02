// --------------------------------------------------------------------------//
/** 
 * AUTHOR: Davide Franco
 * Revised by A. Caminata and S. Marcocci, Sept. 2014 
 */
// --------------------------------------------------------------------------//


//---------------------------------------------------------------------------//

#include "BxLogger.hh"      //Present Bx Class Headers 
#include "BxIO.hh"      //Present Bx Class Headers 

//---------------------------------------------------------------------------//


#include <sstream>
#include <cstdlib>

using namespace std;

static ostringstream devnull;

BxLogger::Severity BxLogger::_minSeverity = BxLogger::debugging;
ostream * BxLogger::_myOstream  = 0 ;
ostream * BxLogger::_myErrstream (&cerr);
ostream * BxLogger::_myNullstream (&devnull);
bool BxLogger::_doPrint = true;
bool BxLogger::_abort=false;

BxLogger::BxLogger(){}

BxLogger::~BxLogger(){}

ostream&  BxLogger::msg(BxLogger::Severity severity, 
		       const char* facility )
{
  _doPrint = true;
  if(!_myOstream) _myOstream = new ostringstream ;
  else ((ostringstream*) _myOstream)->str("");
  
  if(severity >= _minSeverity){
    *_myOstream << toString(severity) << ":" << facility << ":";
  }else{
    _doPrint =false;
    return *_myNullstream ;
  }
/*substituted by creating the new static variable _abort, in order to have the transcription of the error in
 * the log file if severity==fatal. I leave this piece of code here if it turns out the the new one creates problems.
  if(severity == fatal ){
	   *_myOstream  << ::endlog;
    ::abort();
  }*/

  if (severity==fatal)
		  _abort=true;

  return *_myOstream;
}

void  BxLogger::endlog(std::ostream & os){
os<<"";
  if(_doPrint) {
    cout << ((ostringstream*) _myOstream)->str() << endl;
    BxIO::Get()->GetStreamLogFile() << ((ostringstream*) _myOstream)->str() << endl;
  }
  if (_abort==true)
	  ::abort();
}

const char* BxLogger::toString(BxLogger::Severity sever){
  switch (sever) {
  case -2:
    return "Debug";
    break;
  case -1:
    return "Develop";
    break;
  case 0:
    return "Trace";
    break;
  case 1:
    return "Routine";
    break;
  case 2:
    return "Warning";
    break;
  case 3:
    return "Error";
    break;
  case 4:
    return "Fatal";
    break;
  }
  return "Fatal";
}

BxLogger::Severity BxLogger::toEnum(const  std::string& level){
  if(level == "development") return development ;
  if(level == "debugging")   return debugging ;
  if(level == "trace")       return trace;
  if(level == "routine")     return routine;
  if(level == "fatal")       return fatal;
  if(level == "warning")     return warning;
  if(level == "error")       return error;
  return fatal;
}


ostream& endlog(ostream& os){
  BxLogger::endlog( os);
  return os;
}
