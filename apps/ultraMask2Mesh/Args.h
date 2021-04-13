#ifndef ARGS_H
#define ARGS_H

#include <arguments/ArgumentParser.h>
#include "Options.hh"

class Args
{
public:

    /**
     * @brief Args
     * Constructor
     * @param argc
     * Arguments count.
     * @param argv
     * Arguments
     */
    Args(const int argc, const char** argv,  const std::string &help = "");

    /**
     * @brief addArgument
     * Adds a new argument.
     * @param argument
     * A new argument.
     */
    void addArgument(Ultraliser::Argument* argument);

    /**
     * @brief getStringValue
     * Returns the string value that corresponds to the argument.
     * @param argument
     * @return
     */
    std::string getStringValue(Ultraliser::Argument* argument);

    /**
     * @brief getIntegrValue
     * @param argument
     * @return
     */
    int getIntegrValue(Ultraliser::Argument* argument);

    /**
     * @brief getFloatValue
     * @param argument
     * @return
     */
    float getFloatValue(Ultraliser::Argument* argument);

    /**
     * @brief getBoolValue
     * @param argument
     * @return
     */
    bool getBoolValue(Ultraliser::Argument* argument);

    /**
     * @brief parse
     * Parse the command line options.
     */
    void parse();

private:

    /**
     * @brief _parser
     * Command lines arguments parser.
     */
    Ultraliser::ArgumentParser* _parser;

    /**
     * @brief _options
     * Executable options will be encapsulated in this object.
     */
    Options* _options;
};

#endif // ARGS_H
