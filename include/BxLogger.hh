// --------------------------------------------------------------------------//
/** 
 * AUTHOR: Davide Franco
 *Revised by A. Caminata and S. Marcocci, Sept. 2014 
 */
// --------------------------------------------------------------------------//

#ifndef _BXLOGGER_HH
#define _BXLOGGER_HH

#include <iostream>
#include <string>

#define ERRLINE_HACK_1(line)   #line
#define ERRLINE_HACK_2(line)   ERRLINE_HACK_1(line)

#ifdef BxLog
#undef BxLog
#endif
#define BxLog(sev) BxLogger::msg( BxLogger::sev, __FILE__ "(" ERRLINE_HACK_2(__LINE__) ")" )


//---------------------------------------------------------------------------//
/**
 * This class manages the log of the simulation according to the severity
*/
 class BxLogger 
{
public:
/** a severity enum is defined; use only these values
*
*  	fatal:		The message is related to a condition preventing
* 			further execution of the program.  ErrLogger will
*			terminate the program.  Programmers should not call
*			abort or exit themselves.
*
*	error:		A condition exists such that requested result
*			or action can not be produced.  This is serious error.
*
*	warning:	The result is produced, but may not be
* 			what's desired due to an unexpected condition.
*
*	routine:	Nothing known to be wrong with the result;
*			messages that are always produced in normal
*			operation.
*
*	trace:		Messages about the flow of program control
*			and which optional operations took place.
*			(This is the default if nothing is defined)
*  
*	debugging:	Information in addition to the above.
*/
///The severity of the log
  enum Severity {debugging=-2, development=-1, trace=0, routine, warning, error, fatal};

 // members
  static Severity GetSeverity() { return _minSeverity; }
 
/**
 * It writes on the log file the severity level of the message
 * If severity is fatal, it aborts the execution.
 */
 static std::ostream& msg( BxLogger::Severity severity, 
		       const char* facility);//, const char* code );
///It sets the severity
  static void SetSeverity(Severity sever){ _minSeverity = sever;}
 static void endlog(std::ostream &);
  static Severity toEnum(const std::string& );

protected:
 ///default constructor
  BxLogger();

  ///destructor
  ~BxLogger();

  //private  members
private:
  
  static const char* toString(Severity);

  static std::ostream* _myOstream;
  static std::ostream* _myErrstream;
  static std::ostream* _myNullstream;
  static bool _abort;
  static Severity _minSeverity;

  static bool _doPrint;
};
std::ostream& endlog(std::ostream & );

#endif
