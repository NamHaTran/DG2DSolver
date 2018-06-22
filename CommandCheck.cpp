#include "CommandCheck.h"
#include <string>

namespace postProcessKey
{
	std::string exitKey1("Exit"), exitKey2("exit"), exitKey3("EXIT"), exitKey4("End"), exitKey5("end"), exitKey6("END");
	
	/*Function return true if Exit is available*/
	bool checkExit(std::string cmd)
	{
		bool trigger(false);
		if ((cmd.compare(exitKey1) == 0) || (cmd.compare(exitKey2) == 0) || (cmd.compare(exitKey3) == 0) || (cmd.compare(exitKey4) == 0) || (cmd.compare(exitKey5) == 0) || (cmd.compare(exitKey6) == 0))
		{
			trigger = true;
		}
		return trigger;
	}
}

namespace preProcessKey
{
	std::string UnvRdKeyW1("UnvToDG"), UnvRdKeyW2("unvtodg"), UnvRdKeyW3("UNVTODG");
	std::string UnvHelpKeyW1("UnvReaderHelp"), UnvHelpKeyW2("unvreaderhelp"), UnvHelpKeyW3("unvreaderhelp");
	std::string reSubmitKW1("resubmit"), reSubmitKW2("reSubmit");

	/*Function return true if UnvToDG is available*/
	bool checkUnvReader(std::string cmd)
	{
		bool trigger(false);
		if ((cmd.compare(UnvRdKeyW1) == 0) || (cmd.compare(UnvRdKeyW2) == 0) || (cmd.compare(UnvRdKeyW3) == 0))
		{
			trigger = true;
		}
		return trigger;
	}

	/*Function return true if UnvReaderHelp is available*/
	bool checkUnvHelper(std::string cmd)
	{
		bool trigger(false);
		if ((cmd.compare(UnvHelpKeyW1) == 0) || (cmd.compare(UnvHelpKeyW2) == 0) || (cmd.compare(UnvHelpKeyW3) == 0))
		{
			trigger = true;
		}
		return trigger;
	}

	bool reSubmit(std::string cmd)
	{
		bool trigger(false);
		if ((cmd.compare(reSubmitKW1) == 0) || (cmd.compare(reSubmitKW1) == 0))
		{
			trigger = true;
		}
		return trigger;
	}
}

namespace processKey
{
	std::string runningKey1("DG2D"), runningKey2("dg2d");

	/*Function return true if DG2D is available*/
	bool checkDGRun(std::string cmd)
	{
		bool trigger(false);
		if ((cmd.compare(runningKey1) == 0) || (cmd.compare(runningKey2) == 0))
		{
			trigger = true;
		}
		return trigger;
	}
}
