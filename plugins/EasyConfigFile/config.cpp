#include "main.h"
#include "config.h"

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

 string filename;

double double_value;

 int int_value;


void ReadConfig(){
  cout <<"Read config file:"<<endl;
    cout <<"----------------------"<<endl;

    std::ifstream cFile ("config.txt");
    if (cFile.is_open())
    {
        std::string line;
        while(getline(cFile, line)){
            auto delimiterPos = line.find("=");
            auto name = line.substr(0, delimiterPos);
            string value = line.substr(delimiterPos + 1);
            std::cout << name << " " << value << '\n';
       
	    if(name=="file")runname.append(value);
	    if(name=="MyDouble")double_value=stod(value);
	    if(name=="MyInt")int_value=stoi(value);

        }
           
        }
    
    else {
        std::cerr << "Couldn't open config file for reading.\n";
    }
    cout <<"----------------------"<<endl;
  
}

