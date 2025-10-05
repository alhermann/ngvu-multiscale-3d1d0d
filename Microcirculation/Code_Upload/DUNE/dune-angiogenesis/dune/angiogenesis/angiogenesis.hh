#pragma once

#include "includes.hh"
#include "source/RunArtificialNetwork.hh"

class DuneAngiogenesis {
   public:
    DuneAngiogenesis() {}

    virtual ~DuneAngiogenesis() {}

    void run(std::string scenario) {
        std::cout << "===============================================" << std::endl;
        std::cout << "Run Scenario: " << scenario << std::endl;
        std::cout << "===============================================" << std::endl;
        std::cout << " " << std::endl;
        std::cout << " " << std::endl;

        if (scenario == "ArtificialNetwork") {
            RunArtificialNetwork runner;
            runner.run();
        } else {
            std::cout << "===============================================" << std::endl;
            std::cout << "Scenario not implemented!" << std::endl;
            std::cout << "===============================================" << std::endl;
            std::cout << " " << std::endl;
            std::cout << " " << std::endl;
        }
    }
};
