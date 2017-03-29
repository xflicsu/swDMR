/*****************************************************************************
  bedtools.cpp

  bedtools command line interface.  
  Thanks to Heng Li, as this interface is inspired and 
  based upon his samtools interface.

  (c) 2009-2011 - Aaron Quinlan
  Quinlan Laboratory
  Department of Public Health Sciences
  Center for Public Health genomics
  University of Virginia
  aaronquinlan@gmail.com
  
  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "bedtools"

// colors for the term's menu 
#define RESET "\033[m"
#define GREEN "\033[1;32m"
#define BLUE "\033[1;34m"
#define RED "\033[1;31m"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

int intersect_main(int argc, char* argv[]); //
int bedtools_help(void);
int bedtools_faq(void);


int main(int argc, char *argv[])
{
    // make sure the user at least entered a sub_command
    if (argc < 2) return bedtools_help();

    std::string sub_cmd = argv[1];

    // genome arithmetic tools
    if (sub_cmd == "intersect")        return intersect_main(argc-1, argv+1);

    // help
    else if (sub_cmd == "-h" || sub_cmd == "--help" ||
             sub_cmd == "-help")
        return bedtools_help();

    // frequently asked questions
    else if (sub_cmd == "--FAQ" || sub_cmd == "--faq" ||
             sub_cmd == "-FAQ"  || sub_cmd == "-faq")
        return bedtools_faq();

    // verison information
    else if (sub_cmd == "-version" || sub_cmd == "--version")
        cout << "bedtools " << VERSION << endl;

    // verison information
    else if (sub_cmd == "-contact" || sub_cmd == "--contact")
    {
        cout << endl;
        cout << "- For further help, or to report a bug, please " << endl;
        cout << "  email the bedtools mailing list: " << endl;
        cout << "     bedtools-discuss@googlegroups.com" << endl << endl;

        cout << "- Stable releases of bedtools can be found at: " << endl;
        cout << "     http://bedtools.googlecode.com" << endl << endl;

        cout << "- The development repository can be found at: " << endl;
        cout << "     https://github.com/arq5x/bedtools" << endl << endl;
    }
    // unknown
    else {
        // TODO: Implement a Levenstein-based "did you mean???"
        cerr << "error: unrecognized command: " << argv[1] << endl << endl;
        return 1;
    }
    return 0;
}

int bedtools_help(void)
{
    cout  << PROGRAM_NAME  << ": flexible tools for genome arithmetic and DNA sequence analysis.\n";
    cout << "usage:    bedtools <subcommand> [options]" << endl << endl;

    cout  << "The bedtools sub-commands include:" << endl;
    
    cout  << endl;
    cout  << "[ Genome arithmetic ]" << endl;
    cout  << "    intersect     "  << "Find overlapping intervals in various ways.\n";

    cout  << endl;
    cout  << "[ General help ]" << endl;
    cout  << "    --help        "  << "Print this help menu.\n";
    //cout  << "    --faq         "  << "Frequently asked questions.\n";  TODO
    cout  << "    --version     "  << "What version of bedtools are you using?.\n";
    cout  << "    --contact     "  << "Feature requests, bugs, mailing lists, etc.\n";

    cout << "\n";
    return 0;
}


int bedtools_faq(void)
{
    cout << "\n";

    cout << "Q1. How do I see the help for a given command?" << endl;
    cout << "A1. All BEDTools commands have a \"-h\" option. Additionally, some tools " << endl;
    cout << "    will provide the help menu if you just type the command line " << endl;
    cout << "    followed by enter. " << endl;

    cout << "\n";
    return 0;
}
