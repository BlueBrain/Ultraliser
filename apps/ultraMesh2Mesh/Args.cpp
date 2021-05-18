#include "Args.h"
#include <memory>

namespace Ultraliser
{

Args::Args(const int argc, const char **argv, const std::string &help)
{
    // Argument parser
    _parser = new Ultraliser::ArgumentParser(argc, argv, help);
}

Args::~Args()
{
    _parser->~ArgumentParser();
}

void Args::addArgument(Ultraliser::Argument* argument)
{
    _parser->addArgument(argument->getCopy());
}

void Args::parse()
{
    // Show help just in case there is any error before the actual parsing
    _parser->showHelpIfNeeded();

    // Parse the arguments
    _parser->parse();
}

std::string Args::getStringValue(Ultraliser::Argument* argument)
{
    return _parser->getStringValue(argument);
}

int32_t Args::getIntegrValue(Ultraliser::Argument* argument)
{
    return _parser->getIntegrValue(argument);
}

uint32_t Args::getUnsignedIntegrValue(Ultraliser::Argument* argument)
{
    return _parser->getUnsignedIntegrValue(argument);
}

float Args::getFloatValue(Ultraliser::Argument* argument)
{
    return _parser->getFloatValue(argument);
}

bool Args::getBoolValue(Ultraliser::Argument* argument)
{
    return _parser->getBoolValue(argument);
}

}
