#include <iostream>
#include <unistd.h>
#include <iomanip>
#include <sstream>
#include "arg_parse.h"

ArgParse::ArgParse(std::string desc) {
    desc_ = desc;
}

ArgParse::~ArgParse() = default;

std::string ArgParse::get_help(){
    std::stringstream ss;

    ss << "Trans2Genome Help Page"<<std::endl<<std::endl;
    ss << "salmon2genome ";
    for(auto & arg : args_){ // process required arguments
        if(arg.second.required){
            ss<<"-"<<arg.first<<" ";
        }
    }
    bool first = true;
    for(auto & arg : args_){ // now process non-required arguments
        if(!arg.second.required){
            if(first){
                ss<<"[";
            }
            first = false;
            ss<<"-"<<arg.first<<" ";
        }
    }
    if(!first){
        ss<<"]";
    }
    ss<<std::endl;
    ss<<"Arguments:"<<std::endl;
    for(auto & arg : args_){
        ss<<"\t"<<arg.first<<"/--"<<arg.second.name<<"\t"<<arg.second.desc<<std::endl;
    }
    return ss.str();
}

std::string ArgParse::get_param_str() {
    const std::string DELIM = ".";
    std::stringstream ss;

    for (auto a = args_.begin(); a != args_.end(); a++) {
        std::cout<<a->first<<std::endl;

        if (a->second.type != Type::STRING) {
            if (a != args_.begin()) {
                ss << DELIM;
            }
            ss << a->first;

            switch(a->second.type) {
                case Type::FLAG:
                ss << get_flag(a->first);
                break;

                case Type::INT:
                ss << get_int(a->first);
                break;
                
                case Type::DOUBLE:
                ss << std::fixed << std::setprecision(2) << get_double(a->first);
                break;

                default:
                std::cerr << "Error: bad type\n";
                break;
            }
        }
    }
    return ss.str();
}

void ArgParse::parse_args(int argc, char **argv) {

    std::string opt_str;

    for (auto a = args_.begin(); a != args_.end(); a++) {
        opt_str.push_back(a->first);

        if (a->second.type != Type::FLAG) {
            opt_str.push_back(':');
        }
    }

    char o;

    while ( (o = getopt(argc, argv, opt_str.c_str())) != -1 ) {

        if (args_.count(o) == 0) {
            std::cerr << "Error: unrecognized argument " << o << "\n";
            continue;
        }

        
        Arg &a = args_[o];
        a.set = true; // this flag is set in the command line
        switch(a.type) {
            case Type::FLAG:
            *( (bool *) a.value ) = true;
            break;

            case Type::INT:
            *( (int *) a.value ) = atoi(optarg);
            break;
            
            case Type::DOUBLE:
            *( (double *) a.value ) = atof(optarg);
            break;

            case Type::STRING:
            *( (std::string *) a.value ) = std::string(optarg);
            break;
            
            default:
            std::cerr << "Error: bad type\n";
            break;
        }
    }

    // now check if all required arguments have been provided
    for (auto & arg : args_) {
        if(arg.second.required && !arg.second.set){
            std::cerr<<"======================================"<<std::endl<<std::endl<<"Missing argument: "<<arg.second.name<<std::endl<<std::endl<<"======================================"<<std::endl<<std::endl;
            std::cerr<<this->get_help()<<std::endl;
            exit(1);
        }
    }
}

bool ArgParse::add_flag(char c, std::string name, std::string desc="",bool required=false) {

    if (args_.count(c) > 0) {
        return false;
    }

    bool *val_bool = new bool;
    *val_bool = false;

    Arg a = {Type::FLAG, std::move(name), std::move(desc), (void *) val_bool,required,false};
    args_[c] = a;

    return true;
}

bool ArgParse::add_int(char c, std::string name, 
                       int def, std::string desc="",bool required=false) {

    if (args_.count(c) > 0) {
        return false;
    }

    int *val_int = new int;
    *val_int = def;

    Arg a = {Type::INT, std::move(name), std::move(desc), (void *) val_int,required,false};
    args_[c] = a;

    return true;
}

bool ArgParse::add_double(char c, std::string name, 
                          double def, std::string desc="",bool required=false) {

    if (args_.count(c) > 0) {
        return false;
    }

    double *val_double = new double;
    *val_double = def;

    Arg a = {Type::DOUBLE, std::move(name), std::move(desc), (void *) val_double,required,false};
    args_[c] = a;

    return true;

}

bool ArgParse::add_string(char c, std::string name, 
                          std::string def, std::string desc="",bool required=false) {
    
    if (args_.count(c) > 0) {
        return false;
    }

    std::string *val_string = new std::string;
    *val_string = std::move(def);

    Arg a = {Type::STRING, std::move(name), std::move(desc), (void *) val_string,required,false};
    args_[c] = a;

    return true;
}

std::string ArgParse::get_name(char c) {
    return args_[c].name;
}

std::string ArgParse::get_desc(char c) {
    return args_[c].desc;
}

bool ArgParse::get_flag(char c) {
    Arg &a = args_[c];
    return *( (bool *) a.value );
}

int ArgParse::get_int(char c) {
    Arg &a = args_[c];
    return *( (int *) a.value );
}

double ArgParse::get_double(char c) {
    Arg &a = args_[c];
    return *( (double *) a.value );
}

std::string ArgParse::get_string(char c) {
    Arg &a = args_[c];
    return *( (std::string *) a.value );
}

bool ArgParse::is_set(char c){
    Arg &a = args_[c];
    return a.set;
}