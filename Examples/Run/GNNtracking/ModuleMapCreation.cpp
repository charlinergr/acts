// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <memory>
#include "ActsExamples/Framework/Sequencer.hpp"

#include "ActsExamples/GenericDetector/GenericDetector.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"

#include "ActsExamples/Io/Csv/CsvOptionsReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvHitReader.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include "ActsExamples/Io/GraphCreation/GraphCreationOptions.hpp"
#include "ActsExamples/Io/GraphCreation/ModuleMapCreatorWriter.hpp"
#include "ActsExamples/Io/GraphCreation/GraphCreatorWriter.hpp"


//using namespace Acts::UnitLiterals;
using namespace ActsExamples;

int main(int argc, char* argv[]) {
  GenericDetector detector;

  // setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addGeometryOptions(desc);
  Options::addMaterialOptions(desc);
  Options::addInputOptions(desc);
  Options::addGraphCreationOptions(desc);
  Options::addOutputOptions(desc,OutputFormat::Root);

  detector.addOptions(desc);

  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  Sequencer sequencer(Options::readSequencerConfig(vm));

  // Read some standard options
  auto logLevel = Options::readLogLevel(vm);
  auto inputDir = vm["input-dir"].as<std::string>();
  auto outputDir = ensureWritableDirectory(vm["output-dir"].as<std::string>());
  
  // Setup detector geometry
  auto geometry = Geometry::build(vm, detector);
  auto trackingGeometry = geometry.first;
  // Add context decorators
  for (auto cdr : geometry.second) {
    sequencer.addContextDecorator(cdr);
  }

    // Read hits collection
    auto hitReader = Options::readCsvHitReaderConfig(vm);
    //hitReader.outputHits_map = "hits";
    //hitReader.outputHits_vec = "hitsPID";
    sequencer.addReader(
    std::make_shared<CsvHitReader>(hitReader, logLevel));

    //Module Map Creator
    auto minPt = vm["min-pt-cut"].as<float>();
    auto minNHits = vm["min-nhits"].as<long unsigned int>();
    auto rootName = vm["root-filename"].as<std::string>();
    auto giveCutsValues = vm["give-cut-values"].as<bool>();
    // Module map creator script
    ModuleMapCreatorWriter::Config ModuleMapWriterConfig;
    //ModuleMapWriterConfig.trackingGeometry = trackingGeometry;
    ModuleMapWriterConfig.inputDir=inputDir;
    ModuleMapWriterConfig.minPt = minPt;
    ModuleMapWriterConfig.outputDir = outputDir;
    ModuleMapWriterConfig.minNHits = minNHits;
    ModuleMapWriterConfig.rootName = rootName;
    ModuleMapWriterConfig.giveCutsValues = giveCutsValues;
    sequencer.addWriter(
        std::make_shared<ModuleMapCreatorWriter>(ModuleMapWriterConfig, logLevel));
  
  return sequencer.run();
}
