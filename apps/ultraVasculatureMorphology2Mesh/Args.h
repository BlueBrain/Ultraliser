#ifndef ARGS_H
#define ARGS_H

#include <arguments/Arguments.h>
#include <system/System.h>

namespace Ultraliser
{

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
    Args(const int argc, const char** argv, const std::string &help = "");

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
    int32_t getIntegrValue(Ultraliser::Argument* argument);

    /**
     * @brief getUnsignedIntegrValue
     * @param argument
     * @return
     */
    uint32_t getUnsignedIntegrValue(Ultraliser::Argument* argument);

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

}

#endif // ARGS_H
