/*
 *    Part of RMAP software
 *
 *    Copyright (C) 2008 Cold Spring Harbor Laboratory, 
 *                       University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "OptionParser.hpp"

#include <sstream>
#include <cstdlib>
#include <fstream>
#include <cassert>
#include <iomanip>

#include "rmap_utils.hpp"

using std::vector;
using std::string;
using std::endl;

static const size_t MAX_LINE_LENGTH = 72;

enum {
  RMAP_ARG_INT,    RMAP_ARG_UINT,    RMAP_ARG_LONG,
  RMAP_ARG_ULONG,  RMAP_ARG_FLOAT,   RMAP_ARG_DOUBLE,
  RMAP_ARG_STRING, RMAP_ARG_BOOL, RMAP_ARG_CHAR
};

using std::cerr;

void
Option::format_option(const string &argument) {
  std::istringstream ss(argument);
  if ((arg_type == RMAP_ARG_INT && !(ss >> *int_value)) ||
      (arg_type == RMAP_ARG_UINT && !(ss >> *uint_value)) ||
      (arg_type == RMAP_ARG_LONG && !(ss >> *long_value)) ||
      (arg_type == RMAP_ARG_ULONG && !(ss >> *ulong_value)) ||
      (arg_type == RMAP_ARG_FLOAT && !(ss >> *float_value)) ||
      (arg_type == RMAP_ARG_DOUBLE && !(ss >> *double_value)) ||
      (arg_type == RMAP_ARG_CHAR && !(ss >> *char_value)))
    throw RMAPOptionException("Invalid argument [" + argument + 
			      "] to option [" + format_option_name() + "]");
  else if (arg_type == RMAP_ARG_STRING)
    *string_value = argument;
  else if (arg_type == RMAP_ARG_BOOL)
    *bool_value = true;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

Option::Option(const string l_name, const char s_name, const string descr, 
	       const bool reqd, int &val) :
  arg_type(RMAP_ARG_INT), long_name(l_name), short_name(s_name), 
  description(descr), required(reqd), specified(false), int_value(&val)  {}

Option::Option(const string l_name, const char s_name, const string descr,
	       const bool reqd, unsigned int &val) :
  arg_type(RMAP_ARG_UINT), long_name(l_name), short_name(s_name), 
  description(descr), required(reqd), specified(false), uint_value(&val)  {}

Option::Option(const string l_name, const char s_name, const string descr,
	       const bool reqd, long &val) :
  arg_type(RMAP_ARG_LONG), long_name(l_name), short_name(s_name), 
  description(descr), required(reqd), specified(false), long_value(&val) {}

Option::Option(const string l_name, const char s_name, const string descr,
	       const bool reqd, unsigned long &val) :
  arg_type(RMAP_ARG_ULONG), long_name(l_name), short_name(s_name), 
  description(descr), required(reqd), specified(false), ulong_value(&val) {}

Option::Option(const string l_name, const char s_name, const string descr,
	       const bool reqd, float &val) :
  arg_type(RMAP_ARG_FLOAT), long_name(l_name), short_name(s_name), 
  description(descr), required(reqd), specified(false), float_value(&val) {}

Option::Option(const string l_name, const char s_name, const string descr,
	       const bool reqd, double &val) :
  arg_type(RMAP_ARG_DOUBLE), long_name(l_name), short_name(s_name), 
  description(descr), required(reqd), specified(false), double_value(&val) {}

Option::Option(const string l_name, const char s_name, const string descr,
	       const bool reqd, string &val) :
  arg_type(RMAP_ARG_STRING), long_name(l_name), short_name(s_name), 
  description(descr), required(reqd), specified(false), string_value(&val) {}

Option::Option(const string l_name, const char s_name, const string descr, 
	       const bool reqd, bool &val) :
  arg_type(RMAP_ARG_BOOL), long_name(l_name), short_name(s_name), 
  description(descr), required(reqd), specified(false), bool_value(&val) {}

Option::Option(const string l_name, const char s_name, const string descr, 
	       const bool reqd, char &val) :
  arg_type(RMAP_ARG_CHAR), long_name(l_name), short_name(s_name), 
  description(descr), required(reqd), specified(false), char_value(&val) {}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

string
Option::format_option_name() const {
  std::ostringstream ss;
  if (short_name != '\0')
    ss << '-' << short_name << ", -" << long_name;
  else ss << "    -" << long_name;
  return ss.str();
}

string
Option::format_option_description(const size_t offset) const {
  std::ostringstream ss;
  if (!description.empty()) {
    vector<string> parts;
    rmap::split_whitespace(description, parts);
    
    size_t line_len = 0;
    for (size_t i = 0; i < parts.size(); ++i) {
      if (offset + line_len + parts[i].size() >= MAX_LINE_LENGTH && i > 0) {
	line_len = 0;
	ss << endl;
      }
      if (i > 0 && line_len == 0)
	ss << string(offset, ' ');
      ss << parts[i] << " ";
      line_len += parts[i].size();
    }
  }
  return ss.str();
}

bool
Option::option_match(const string &other) {
  return (long_name == other || 
	  (other.length() > 1 && other[0] == '-' && 
	   (other.substr(1) == long_name || 
	    (other[1] == short_name && other.length() == 2))));
}

bool
Option::parse(vector<string> &command_line) {
  static const string dummy("");
  if (!command_line.empty())
    for (size_t i = 0; i < command_line.size();)
      if (option_match(command_line[i])) {
	if (i < command_line.size() - 1) 
	  format_option(command_line[i + 1]);
	else format_option(dummy);
	specified = true;
	command_line.erase(command_line.begin() + i);
	if (arg_type != RMAP_ARG_BOOL)
	  command_line.erase(command_line.begin() + i);
      }
      else ++i;
  return (specified || !required);
}

void
Option::parse_config_file(vector<string> &options) {
  for (size_t i = 0; i < options.size(); ++i) {
    vector<string> opt_val = rmap::split(options[i], "=");
    if (option_match(opt_val.front())) {
      format_option(opt_val.back());
      options.erase(options.begin() + i);
      specified = true;
    }
  }
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

void
OptionParser::add_opt(const string l_name, const char s_name, const string descr, 
		      const bool reqd, int &val) {
  options.push_back(Option(l_name, s_name, descr, reqd, val));
}

void
OptionParser::add_opt(const string l_name, const char s_name, const string descr,
		      const bool reqd, unsigned &val) {
  options.push_back(Option(l_name, s_name, descr, reqd, val));
}

void
OptionParser::add_opt(const string l_name, const char s_name, const string descr,
		      const bool reqd, long &val)  {
  options.push_back(Option(l_name, s_name, descr, reqd, val));
}

void 
OptionParser::add_opt(const string l_name, const char s_name, const string descr,
		      const bool reqd, unsigned long &val) {
  options.push_back(Option(l_name, s_name, descr, reqd, val));
}

void 
OptionParser::add_opt(const string l_name, const char s_name, const string descr,
		      const bool reqd, float &val) {
  options.push_back(Option(l_name, s_name, descr, reqd, val));
}

void 
OptionParser::add_opt(const string l_name, const char s_name, const string descr,
		      const bool reqd, double &val) {
  options.push_back(Option(l_name, s_name, descr, reqd, val));
}

void 
OptionParser::add_opt(const string l_name, const char s_name, const string descr, 
		      const bool reqd, string &val) {
  options.push_back(Option(l_name, s_name, descr, reqd, val));
}

void 
OptionParser::add_opt(const string l_name, const char s_name, const string descr,
		      const bool reqd, bool &val) {
  options.push_back(Option(l_name, s_name, descr, reqd, val));
}

void 
OptionParser::add_opt(const string l_name, const char s_name, const string descr,
		      const bool reqd, char &val) {
  options.push_back(Option(l_name, s_name, descr, reqd, val));
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#include <cstring>

static void
read_config_file(const string &config_filename, vector<string> &config_file_options) {
  static const char COMMENT_CHARACTER = '#';
  static const size_t INPUT_BUFFER_SIZE = 1000;
  config_file_options.clear();
  
  std::ifstream in(config_filename.c_str());
  if (!in)
    throw RMAPOptionException("cannot open config file " + config_filename);
  
  size_t line_number = 0;
  while (!in.eof()) {
    ++line_number;
    char buffer[INPUT_BUFFER_SIZE + 1];
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() == static_cast<int>(INPUT_BUFFER_SIZE))
      throw RMAPOptionException("Line in " + config_filename + 
				  "\nexceeds max length: " +
				  toa(INPUT_BUFFER_SIZE));
    if (in.gcount() == 0)
      throw RMAPOptionException("Problem reading file \"" + config_filename +
				  "\" (empty?)");
    const size_t final_char = in.gcount() - 1;
    if (buffer[final_char] == '\r') buffer[final_char] = '\0';
    if (buffer[0] != '\0' && buffer[0] != COMMENT_CHARACTER) {
      string line(buffer);
      line = rmap::strip(line);
      if (line.length() == 0 || line.find('=') == string::npos)
	throw RMAPOptionException("Line " + toa(line_number) +
				    " malformed in config file " + config_filename);
      config_file_options.push_back(line);
    }
    in.peek();
  }
}

void
OptionParser::parse(const int argc, const char **argv, vector<string> &arguments,
		    string config_filename) {
  // The "2" below corresponds to the "about" and "help" options
  assert(options.size() >=  2);

  if (!config_filename.empty()) {
    vector<string> config_file_options;
    read_config_file(config_filename, config_file_options);
    for (size_t i = 0; i < options.size(); ++i)
      options[i].parse_config_file(config_file_options);
  }
  arguments.clear();

  // The '1' and '+ 1' below is to skip over the program name
  assert(argc >= 1);
  copy(argv + 1, argv + argc, back_inserter(arguments));
  
  for (size_t i = 0; i < options.size(); ++i)
    if (!options[i].parse(arguments) && first_missing_option_name.empty())
      first_missing_option_name = options[i].format_option_name();
}

OptionParser::OptionParser(const string nm, const string descr,
			   string noflag_msg) :
  prog_name(nm), prog_descr(descr), noflag_message(noflag_msg),
  help_request(false), about_request(false) {
  add_opt("help", '?', "print this help message", false, help_request);
  add_opt("about", '\0', "print about message", false, about_request);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////// 
////// FOR PRINTING MESSAGES
//////

string 
OptionParser::help_message() const {
  // corresponds to the two spaces before and 
  static const char *SPACE_BEFORE_SHORT = "  ";
  static const char *SPACE_BTWN_SHRT_LNG = "  ";
  static const size_t TOTAL_ADDED_SPACE = 4;
  
  vector<string> option_names;
  size_t max_name_len = 0;
  for(size_t i = 0; i < options.size(); ++i) {
    option_names.push_back(options[i].format_option_name());
    max_name_len = std::max(max_name_len, option_names.back().length());
  }
  
  std::ostringstream ss;
  ss << "Usage: " << prog_name << " [OPTIONS]";
  if (!noflag_message.empty())
    ss << " " << noflag_message;
  ss << endl << endl;
  
  if (options.size() > 2) {
    ss << "Options:" << endl;
    for (size_t i = 2; i < options.size(); ++i)
      ss << SPACE_BEFORE_SHORT << std::left << std::setw(max_name_len) 
	 << option_names[i] << SPACE_BTWN_SHRT_LNG
	 << options[i].format_option_description(max_name_len + 
						 TOTAL_ADDED_SPACE) << endl;
  }
  
  ss << endl << "Help options:" << endl;
  for (size_t i = 0; i < std::min(2ul, options.size()); ++i)
    ss << SPACE_BEFORE_SHORT << std::left << std::setw(max_name_len) 
       << option_names[i] << SPACE_BTWN_SHRT_LNG
       << options[i].format_option_description(max_name_len + 
					       TOTAL_ADDED_SPACE) << endl;
  return ss.str();
}

string 
OptionParser::about_message() const {
  static const char *PROGRAM_NAME_TAG =  "PROGRAM: ";
  
  vector<string> parts;
  rmap::split_whitespace(prog_descr, parts);

  size_t line_len = 0;
  std::ostringstream ss;
  ss << PROGRAM_NAME_TAG << prog_name << endl;
  for (size_t i = 0; i < parts.size(); ++i) {
    if (line_len + parts[i].size() >= MAX_LINE_LENGTH && i > 0) {
      line_len = 0;
      ss << endl;
    }
    ss << parts[i] << " ";
    line_len += parts[i].size();
  }
  return ss.str();
}

string
OptionParser::option_missing_message() const {
  std::ostringstream ss;
  ss << "required argument missing: [" << first_missing_option_name << "]";
  return ss.str();
}
