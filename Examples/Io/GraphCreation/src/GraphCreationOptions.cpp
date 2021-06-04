// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/GraphCreation/GraphCreationOptions.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include <string>
#include <boost/program_options.hpp>

void ActsExamples::Options::addGraphCreationOptions(
    ActsExamples::Options::Description& desc) {
  using boost::program_options::value;

  auto opt = desc.add_options();
  opt("min-pt-cut",value<float>()->default_value(0.),"Min cut on the generated particles");
  opt("min-nhits",value<long unsigned int>()->default_value(2),"Min cut of number of hits of generated particles");
  opt("root-filename",value<std::string>(), "Name of the root file");
  opt("give-cut-values",value<bool>()->default_value("False"), "Give cuts values");
  opt("input-module-map",value<std::string>(), "Name of the root file");
  opt("give-true-graph",value<bool>()->default_value("False"), "Give true graph");
  opt("save-graph-on-disk",value<bool>()->default_value("False"), "Save graph on disk");
}

ActsExamples::ModuleMapCreatorWriter::Config
ActsExamples::Options::readModuleMapWriterConfig(
    const ActsExamples::Options::Variables& variables) {
  auto minPt = variables["min-pt-cut"].template as<float>();
  auto minNHits = variables["min-nhits"].template as<long unsigned int>();
  auto rootName = variables["root-filename"].template as<std::string>();
  auto giveCutsValues = variables["give-cut-values"].template as<bool>();

  ModuleMapCreatorWriter::Config cfg;
  cfg.minPt = minPt;
  cfg.minNHits = minNHits;
  cfg.rootName = rootName;
  cfg.giveCutsValues = giveCutsValues;
  
  return cfg;
}

ActsExamples::GraphCreatorWriter::Config
ActsExamples::Options::readGraphCreatorConfig(
    const ActsExamples::Options::Variables& variables) {
  auto inputModuleMap = variables["input-module-map"].template as<std::string>();
  auto minPt = variables["min-pt-cut"].template as<float>();
  auto minNHits = variables["min-nhits"].template as<long unsigned int>();
  auto trueGraph = variables["give-true-graph"].template as<bool>();
  auto saveGraph = variables["save-graph-on-disk"].template as<bool>();

  GraphCreatorWriter::Config cfg;
  cfg.inputModuleMap = inputModuleMap;
  cfg.trueGraph = trueGraph;
  cfg.saveGraph = saveGraph;
  cfg.minPt = minPt;
  cfg.minNHits = minNHits;
  
  return cfg;
}