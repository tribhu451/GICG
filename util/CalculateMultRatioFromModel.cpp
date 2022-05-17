#include <iostream>
#include <cmath>
#include <fstream>
#include<string>
#include<sstream>

using namespace std;
using std::cout;
using std::cin;
using std::fstream ; 
//using std::sstream ; 
using std::endl;
using std::string;

#define PI 3.141592653

int main()
{

  double RefMult ; 
  double RefMult_err ; 
  cout << "Enter Reference Multiplicity." << endl ;  
  cin  >> RefMult ;
  cout << "Enter Reference Multiplicity error." << endl ;  
  cin  >> RefMult_err ;
 
 
  string input_file_name = "Energy.txt" ; // this is the output file obtained after running the code centrality_determination.C
  string output_file_name = "EnergyRatio.txt" ; 


  istringstream* iss ;
  char   buff[1000] ;
  double dummy0,CentLow,CentHigh,CentMid,DnchDeta,DnchDeta_err ;

  // open the file
  std::fstream InFile;
  InFile.open(input_file_name.c_str(),std::ios::in);
  if(!InFile){
    cout<<"Input file not found."<<endl;
    exit(1);
  } 

   std::fstream OutFile;
   OutFile.open(output_file_name.c_str(),std::ios::out);
   OutFile << "#CentralityLow	CentralityHigh	centrlityMidValue	dNch/deta|eta=0		MultiplicityRatio   MultiplicityRatio_error" << endl ;

   // read the first line  
   InFile.getline(buff,1000);

   for(int i=0; i<11; i++){
    InFile.getline(buff,1000);
    iss = new istringstream(buff);
    *iss >> CentLow >> CentHigh >> CentMid >> dummy0 >> dummy0 >> dummy0 >> dummy0 >> DnchDeta >> DnchDeta_err >> dummy0 ;
    cout << CentLow  << "  " << CentHigh << "  " << CentMid << "  " << DnchDeta << "  "
                << "  " << DnchDeta / RefMult << "  " 
                   << sqrt( pow(DnchDeta,2) / pow(RefMult,4) * pow(RefMult_err,2) + 1.0 / pow(RefMult,2) * pow(DnchDeta_err,2) )   << endl ; 
    OutFile << CentLow  << "  " << CentHigh << "  " << "  " << CentMid << "  " << DnchDeta << "  "
                << "  " << DnchDeta / RefMult << "  "
                  << sqrt( pow(DnchDeta,2) / pow(RefMult,4) * pow(RefMult_err,2) + 1.0 / pow(RefMult,2) * pow(DnchDeta_err,2) )  << endl ; 

    delete iss;  
   }

   OutFile.close();
   InFile.close();
   return 0 ; 

}



