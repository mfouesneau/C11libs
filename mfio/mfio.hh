/**
 * MFIO - Trying to make ascii IO simpler.
 * =======================================
 * 
 * `AsciiReader` and `AsciiWriter` make very simple proxies to reading and writing
 * Ascii data. They both also handle `stdin` and `stdout`.
 * 
 * @author M. Fouesneau
 * 
 * example
 * -------
 * 
 * .. code:: c
 * 
 *          // open data file 
 *          AsciiReader rd("file.csv");
 *          AsciiWriter wt("-");        // will write on stdout
 * 
 *          // set delimiter and comment char;
 *          rd.set_delimiter(",");
 *          rd.set_comment("#");
 * 
 *          // wt.set_delimiter(" "); // space is default value
 * 
 *          // select columns to use
 *          rd.set_columns({0, 1, 3});
 * 
 *          // define data field values and types
 *          double ra;
 *          double dec;
 *          double gmag;
 * 
 *          // while rd can read lines from the input file
 *          while(++rd){
 *                 rd.get_next_value<double>(ra);
 *                 rd.get_next_value<double>(dec);
 *                 rd.get_next_value<double>(gmag);
 *                 wt << {ra, dec, gmag};  // write a vector in one command
 *                 ++ wt;                  // new line
 *          }
 * 
 * Complex string formatting exports can be done through `io::string_format`, which
 * gives the same API as `printf`.
 * 
 * 
 * AsciiReader
 * -----------
 * 
 * Reads in ascii data.
 * 
 * Notes
 * ~~~~~
 * 
 * * `count_lines` is based on GNU wc implementation for speed. It can skip
 *   commented lines and be used to prepare memory arrays. It also caches the
 *   results after the first call.
 * 
 * * reading a new line (`operator ++`) will always skip commented lines.
 *  
 * * parsing values through `get_next_value` will ignore non-used columns (if set).
 */
#pragma once

#include <stdlib.h>
#include <iostream>
#include <sstream>      // std::stringstream
#include <fstream>      // fstream 
#include <vector>       // std::vector
#include <algorithm>
#include <string>
#include <cstdio>
#include <cstdarg>     // for ... as arguments


/**
 * Split a string into a vector given a delimiter.
 *
 * note that consecutive delimiters are processed as one.
 *
 * @param line string to split
 * @param delimiter delimiter of the tokens
 *
 * @return vector of tokens
 */
std::vector<std::string> split_string(const std::string& line, const char delimiter){

    std::vector<std::string> tokens;
    std::string token;

    std::for_each(line.begin(), line.end(), [&](char c){
        if (c != delimiter){
            token += c;
        } else {
            if (token.length()) tokens.push_back(token);
            token.clear();
        }
        });
    if (token.length()) tokens.push_back(token);

    return tokens;
}


namespace io {

    /**
     * Parse a string to a numeric format.
     *
     * @param str string to parse
     *
     * @return value of type T
     */
    template <typename T>
    T parse(const std::string& str){
        T val;
        std::istringstream(str) >> val;
        return val;
    }

    /**
     * Generate a formatted string a la fprintf.
     *
     * @param fmt format to use
     * @param args relevant information 
     *
     * @return string
     */
    inline std::string string_format(const char* fmt, ...){
        int size = 512;
        char* buffer = 0;
        buffer = new char[size];
        va_list vl;
        va_start(vl, fmt);
        int nsize = vsnprintf(buffer, size, fmt, vl);
        if(size<=nsize){ //fail delete buffer and try again
            delete[] buffer;
            buffer = 0;
            buffer = new char[nsize+1]; //+1 for /0
            nsize = vsnprintf(buffer, size, fmt, vl);
        }
        std::string ret(buffer);
        va_end(vl);
        delete[] buffer;
        return ret;
    }


    /**
     * Returns a human readable string from a number of bytes.
     *
     * @param bytes number of bytes
     *
     * @return the corresponding string
     * */
    std::string format_bytes(long bytes)
    {
        std::vector<std::string> suffix = { "B", "KB", "MB", "GB", "TB" };
        size_t i;
        double dblSByte = bytes;
        for (i = 0; (i < suffix.size()) && (bytes >= 1024); ++i, bytes /= 1024){
            dblSByte = bytes / 1024.0;
        }

        return string_format("%.2f ", dblSByte) + suffix[i];
    }


    /**
    * Counts the number of lines on the input buffer
    *
    * @param fname        input file to read from
    * @param buffer_size  the size of the buffer (too large/small is not good.
    *                     16K = 16384 is what GNU wc command uses
    *                     64K is the pipe buffer in Linux systems
    * @param comment      Do not count line if starts with this character
    *
    * @return nlines      the number of lines found in stream
    */
    size_t count_lines(const std::string& fname, const std::string& comment="#", 
            const int buffer_size = 8 * 1024 * 1024) {

        std::FILE* stream = std::fopen(fname.c_str(), "r");
        std::vector<char> buffer(buffer_size);
        size_t size;
        size_t line_count = 0;
        size_t comment_count = 0;
        char prev = '\n';
        bool comment_line = false;
        while ((size = std::fread(buffer.data(), sizeof(char), buffer_size, stream)) > 0) {
        for (size_t i=0; i < size; ++i){
            char* ch = buffer.data() + i;
            if ((*ch == comment[0]) && (prev == '\n')){
                comment_line = true;
            }
            if (*ch == '\n') {
                if (!comment_line) { 
                    line_count ++;
                } 
                else { 
                    comment_line = false;
                }
            }
            prev = *ch;
        }
        // tmp = std::count_if(buffer.begin(), 
        //                    buffer.begin() + size, 
        //                    [](char ch) { return ch == '\n'; });
        }
        std::fclose(stream);
        return line_count - comment_count;
    }

    /**
    * Counts the number of lines on the input buffer.
    *
    * note that this version is optimize if no comment char was provided.
    *
    * @param fname        input file to read from
    * @param buffer_size  the size of the buffer (too large/small is not good.
    *                     16K = 16384 is what GNU wc command uses
    *                     64K is the pipe buffer in Linux systems
    *
    * @return nlines      the number of lines found in stream
    */
    size_t count_lines(const std::string fname, const int buffer_size = 8 * 1024 * 1027){
        std::FILE* stream = std::fopen(fname.c_str(), "r");
        std::vector<char> buffer(buffer_size);
        size_t size;
        size_t line_count = 0;
        while ((size = std::fread(buffer.data(), sizeof(char), buffer_size, stream)) > 0) {
            line_count += std::count_if(buffer.begin(), 
                    buffer.begin() + size, 
                    [](char ch) { return ch == '\n'; });
        }
        return line_count;
    }


    /**
     * Buffer the next line from a stream ignoring comments.
     *
     * @param stream  stream to read from
     * @param line    line to write into
     * @param comment comment character
     *
     * @return true if anything was loaded.
     */
    bool get_next_line(std::ifstream& stream, std::string& line, 
            const std::string comment="#"){

        bool val = false;
        line.clear();
        while(std::getline(stream, line)){
            // skip line if comment line
            if (comment.find_first_of(line) == 0){
                continue;
            }
            val = true;
            break;
        }
        return val;
    }


    /**
     * Return the first non comment line from the file
     *
     * @param fname   file to open and read from
     * @param line    line to write into
     * @param comment comment character
     *
     * @return vector of values, ie, field names
     */
    std::vector<std::string> get_header(const std::string fname, 
            const std::string comment="#",
            const std::string delimiter=","){

        std::ifstream stream(fname);
        std::string line;
        get_next_line(stream, line, comment);
        stream.close();
        return split_string(line, delimiter[0]);
    }


    /**
     * Convenient class to use individual functions
     * Works also on stdin
     */
    class AsciiReader {

    public:
    	AsciiReader (const std::string& fname);
	    size_t count_lines(const int buffer_size = 8 * 1024 * 1024);
    	bool is_column_used(size_t col);
    	bool operator++();
	
	    // getters
    	std::vector<std::string> get_header();
    	std::string get_comment();
    	std::string get_delimiter();
    	std::vector<size_t> get_columns();
        std::size_t get_number_of_fields();
	    template <typename T> bool get_next_value(T& val);

	    // setters
    	void set_comment(const std::string& comment);
    	void set_delimiter(const std::string& delimiter);
    	void set_columns(const std::vector<size_t>& usecols);

    private:
    	/* data */
        std::string fname;                         /**< filename */
        std::string comment = "#";                 /**< comment character */
        std::string delimiter = " ";               /**< delimiter character */
        std::vector<size_t> usecols = {};          /**< column to use only */

    	std::shared_ptr<std::istream> is;          /**< pointer to stream */
        std::string current_line = "";             /**< current line from is */
        std::vector<std::string> current_values;   /**< current line from is */
        size_t current_column = 0;                 /**< current column position */
        size_t nlines = 0;                         /**< cache number of lines */
        size_t nfields = 0;                        /**< cache number of fields */
        std::vector<std::string> header;           /**< cache header */


        void reset_stream(const std::string& fname);
    };

    /**
     * Simple constructor
     */
    AsciiReader::AsciiReader(const std::string& fname){
	    // std::cout << "AsciiReader: opening " << fname << std::endl;
	    this->reset_stream(fname);
    }

    /** 
     * Set columns to use.
     *
     * @param usecols indexes of columns to use
     *
     * @throws runtime_error if columns are outside possible range
     */
    void AsciiReader::set_columns(const std::vector<size_t>& usecols){
	    this->usecols = std::vector<size_t>(usecols);
	    std::sort(this->usecols.begin(), this->usecols.end());
	    if (usecols.back() > this->get_header().size())
		    throw std::runtime_error("Column index above maximum number of fields");
    }

    /** 
     * get columns to use.
     *
     * @return usecols indexes of columns to use
     */
    std::vector<size_t> AsciiReader::get_columns(){
	    return this->usecols;
    }


    /**
     * returns how many columns will be read
     *
     * @return this->get_header().size() or number of selected columns
     */
    size_t AsciiReader::get_number_of_fields(){
        if (this->nfields > 0)
            return this->nfields;

        size_t n = this->usecols.size();
        if (n)
            this->nfields = n;
        else
            this->nfields = this->get_header().size();

        return this->nfields;
    }

    /** 
     * Set comment string
     *
     * @param comment string that starts a comment line
     */
    void AsciiReader::set_comment(const std::string& comment){
	    this->nlines = -1; // reset cache
	    this->comment = comment;
    }

    /** 
     * Get comment string
     *
     * @return comment string that starts a comment line
     */
    std::string AsciiReader::get_comment(){
	    return this->comment;
    }


    /** 
     * Set delimiter string
     *
     * @param delimiter string that starts a delimiter line
     */
    void AsciiReader::set_delimiter(const std::string& delimiter){
	    this->nlines = 0; // reset cache
	    this->delimiter = delimiter;
    }

    /** 
     * Get delimiter string
     *
     * @return delimiter string that starts a delimiter line
     */
    std::string AsciiReader::get_delimiter(){
	    return this->delimiter;
    }


    /**
     * open the stream to the file or stdin.
     *
     * @param fname  filename (empty or "-" means stdin)
     *
     * @throws runtime_error if the file cannot be opened.
     */
    void AsciiReader::reset_stream(const std::string& fname){
	    if ((fname.size() == 0) or (fname.compare("-") == 0)){
            // use stdin
            this->is.reset(&std::cin, [](...){}); // do nothing for delete
            this->fname = "stdin";
	    } else {
            this->fname = fname;
            this->is.reset(new std::ifstream(fname));
            if (!*this->is){
                 throw std::runtime_error("\nERROR: Cannot open file " + fname);
            }
	    }
    }


    /**
     * Counts the number of non-commented lines.
     *
     * @param buffer_size  the size of the buffer (too large/small is not good.
     *                     16K = 16384 is what GNU wc command uses
     *                     64K is the pipe buffer in Linux systems
     *
     * @return nlines      the number of lines found in stream
     */
    size_t AsciiReader::count_lines(const int buffer_size){

	    // no counting on stdin or stream
	    if (this->fname == "stdin"){
		    return -1;
	    }

	    if (this->nlines > 0){
		    return this->nlines;
	    }
	    if (this->comment.size()) {
		    this->nlines = io::count_lines(this->fname, this->comment, buffer_size);
	    } else {
		    // use the fast version
		    this->nlines = io::count_lines(this->fname, buffer_size);
	    }
	    return this->nlines;
    }

    /**
     * Load next non-comment line (if comment is defined)
     *
     * @return true is something was loaded
     */
    bool AsciiReader::operator++(){
        bool val = false;
        this->current_line.clear();
	    this->current_column = 0;
	    
	    if (this->comment.size()){
		    while(std::getline(*this->is, this->current_line)){
                // skip line if comment line
                if (this->comment.find_first_of(this->current_line) == 0){
                    continue;
                }
                val = true;
                break;
		    }
	    } else {
		    while(std::getline(*this->is, this->current_line)){
			    val = true;
			    break;
		    }

	    }
        if (current_line.size()){
            this->current_values = split_string(this->current_line, this->delimiter[0]);
        }
        return val;
    }


    /**
     * Return the first non comment line from the file
     *
     * @param fname   file to open and read from
     * @param line    line to write into
     * @param comment comment character
     *
     * @return vector of values, ie, field names
     */
    std::vector<std::string> AsciiReader::get_header(){ 
     
        if (this->fname == "stdin"){
            return std::vector<std::string>();  // empty
        }
        // cache
        if (this->header.size()){
            return this->header;
        }

        std::ifstream stream(this->fname);
        std::string line;
        io::get_next_line(stream, line, comment);
        stream.close();
        this->header = split_string(line, delimiter[0]);
        return this->header;
    }


    /**
     * Check if a column index is used.
     *
     * @param col  column index to check
     *
     * @return true if the column is used
     */
    bool AsciiReader::is_column_used(size_t col){
	// if there is no pre-defined column set then use them all
	if (this->usecols.size() == 0){
		return true;
	}
	// scan through the column set
	// if col is in, use it
	// columns are sorted so that if the value is above col stop looking
	// too.
    if (this->usecols.back() < col){
        return false;
    }
    for (size_t key=0; key < this->usecols.size(); ++key){
        if (this->usecols[key] == col){
            return true;
                } else if (this->usecols[key] > col){
            return false;
        }
    }
	return false;
    }


    /**
     * Read the next value while skipping not used columns
     *
     * @param val value into which store the data
     *
     * @return true is something was valid
     */
    template <typename T>
    bool AsciiReader::get_next_value(T& val){
	
        // Current column position is used
        if (this->is_column_used(this->current_column)){
            if (this->current_column >= this->current_values.size()) 
                return false;
            std::string sfield = this->current_values[this->current_column];
            if (sfield.size()){
                val = parse<T>(sfield);
                ++ this->current_column;
                return true;
            }
        } 
        ++ this->current_column;
        return this->get_next_value<T>(val);
    }

    /**
     * A simple class to write ascii outputs on disk (or stdout)
     */
    class AsciiWriter {

    public:
    	AsciiWriter (const std::string& fname);
	
	    // getters
    	std::string get_comment();
    	std::string get_delimiter();

	    // setters
    	void set_comment(const std::string& comment);
    	void set_delimiter(const std::string& delimiter);

        template<typename T> AsciiWriter& operator<< (const T& val);
        template<typename T> AsciiWriter& operator<< (const std::vector<T>& val);
        void operator++();

    private:
    	/* data */
        std::string fname;                         /**< filename */
        std::string comment = "#";                 /**< comment character */
        std::string delimiter = " ";               /**< delimiter character */

    	std::shared_ptr<std::ostream> os;          /**< pointer to stream */
        size_t current_line = 0;                   /**< cache number of lines */
        size_t current_col = 0;                    /**< count number of columns */

        void reset_stream(const std::string& fname);
    };


    /** Constructor
     * 
     * @param fname file name
     */
    AsciiWriter::AsciiWriter(const std::string& fname) {
        this->reset_stream(fname);
    }

    /** Set the stream
     *
     * @param fname file name
     *
     */
    void AsciiWriter::reset_stream(const std::string& fname){
	    if ((fname.size() == 0) or (fname.compare("-") == 0)){
            // use stdout
            this->os.reset(&std::cout, [](...){}); // do nothing for delete
            this->fname = "stdout";
	    } else {
            this->fname = fname;
            this->os.reset(new std::ofstream(fname));
            if (!*this->os){
                 throw std::runtime_error("\nERROR: Cannot open file " + fname);
            }
	    }
    }

    /** 
     * Set comment string
     *
     * @param comment string that starts a comment line
     */
    void AsciiWriter::set_comment(const std::string& comment){
	    this->comment = comment;
    }

    /** 
     * Get comment string
     *
     * @return comment string that starts a comment line
     */
    std::string AsciiWriter::get_comment(){
	    return this->comment;
    }


    /** 
     * Set delimiter string
     *
     * @param delimiter string that starts a delimiter line
     */
    void AsciiWriter::set_delimiter(const std::string& delimiter){
	    this->delimiter = delimiter;
    }

    /** 
     * Get delimiter string
     *
     * @return delimiter string that starts a delimiter line
     */
    std::string AsciiWriter::get_delimiter(){
        return this->delimiter;
    }

    /**
     * Add value to the current output. Adds the delimiter if needed
     *
     * @param val value to add
     *
     * @returns itself for further inputs
     */
    template<typename T>
    AsciiWriter& AsciiWriter::operator<< (const T& val){
        std::ostream *s = this->os.get();
        if (this->current_col > 0){
            *s << this->delimiter;
        }
        *s << val;
        ++(this->current_col);
        return *this;
    }

    /**
     * Add multiple values to the current output. Adds the delimiter when needed
     *
     * @param val vector of values to add
     *
     * @returns itself for further inputs
     */
    template<typename T>
    AsciiWriter& AsciiWriter::operator<< (const std::vector<T>& val){
        std::ostream *s = this->os.get();
        for (T v: val){
            if (this->current_col > 0){
                *s << this->delimiter;
            } 
            *s << v;
            ++(this->current_col);
        }
        return *this;
    }

    /**
     * Set new line in the output
     */
    void AsciiWriter::operator++(){
        std::ostream *s = this->os.get();
        *s << std::endl;
        (this->current_line)++;
    }
}

// vim: expandtab:ts=4:softtabstop=4:shiftwidth=4
