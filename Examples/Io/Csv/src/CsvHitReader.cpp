// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvHitReader.hpp"
#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Utilities/Range.hpp"
#include <map>
#include <boost/bimap.hpp>
#include <dfe/dfe_io_dsv.hpp>
#include "CsvOutputData.hpp"
#include <boost/unordered_map.hpp>
#include "ActsExamples/EventData/SimParticle.hpp"

ActsExamples::CsvHitReader::CsvHitReader(
    const ActsExamples::CsvHitReader::Config& cfg,
    Acts::Logging::Level lvl)
    : m_cfg(cfg),
      // TODO check that all files (hits,cells,truth) exists
      m_eventsRange(determineEventFilesRange(cfg.inputDir, "hits.csv")),
      m_logger(Acts::getDefaultLogger("CsvHitReader", lvl)) {
}

std::string ActsExamples::CsvHitReader::CsvHitReader::name()
    const {
  return "CsvHitReader";
}

std::pair<size_t, size_t> ActsExamples::CsvHitReader::availableEvents()
    const {
  return m_eventsRange;
}

namespace {
    struct CompareHitId {
    // support transparent comparision between identifiers and full objects
    using is_transparent = void;
    template <typename T>
    constexpr bool operator()(const T& left, const T& right) const {
        return left.hit_id < right.hit_id;
    }
    template <typename T>
    constexpr bool operator()(uint64_t left_id, const T& right) const {
        return left_id < right.hit_id;
    }
    template <typename T>
    constexpr bool operator()(const T& left, uint64_t right_id) const {
        return left.hit_id < right_id;
    }
    };

    struct CompareParticleId {
    // support transparent comparision between identifiers and full objects
    using is_transparent = void;
    template <typename T>
    constexpr bool operator()(const T& left, const T& right) const {
        return left.particle_id < right.particle_id;
    }
    template <typename T>
    constexpr bool operator()(uint64_t left_pid, const T& right) const {
        return left_pid < right.particle_id;
    }
    template <typename T>
    constexpr bool operator()(const T& left, uint64_t right_pid) const {
        return left.particle_id < right_pid;
    }
    };

    template <typename Data>
    inline std::vector<Data> readEverything(
        const std::string& inputDir, const std::string& filename,
        const std::vector<std::string>& optionalColumns, size_t event) {
        std::string path = ActsExamples::perEventFilepath(inputDir, filename, event);
        dfe::NamedTupleCsvReader<Data> reader(path, optionalColumns);

        std::vector<Data> everything;
        Data one;
        while (reader.read(one)) {
            everything.push_back(one);
        }
        return everything;
    }

    //template <typename Data>
    inline boost::unordered_map<uint64_t,std::vector<uint64_t>> readEverythingStruct(
        const std::string& inputDir, const std::string& filename,
        const std::vector<std::string>& optionalColumns, size_t event) {
        std::string path = ActsExamples::perEventFilepath(inputDir, filename, event);
        dfe::NamedTupleCsvReader<ActsExamples::TruthHitData> reader(path, optionalColumns);
        boost::unordered_map<uint64_t,std::vector<uint64_t>> map;
        ActsExamples::TruthHitData one;
        while (reader.read(one)) {
            //std::cout<<one.particle_id<<" "<<one.hit_id<<std::endl;
            uint64_t pid = one.particle_id;
            uint64_t hit = one.hit_id;
            map[pid].push_back(hit);
            //map.insert(pid,hit);
        }
        return map;
    }

    std::vector<ActsExamples::HitData> readHits(const std::string& inputDir, size_t event) {
        auto hits = readEverything<ActsExamples::HitData>(inputDir, "hits.csv", {}, event);
        return hits;
    }

    std::vector<ActsExamples::TruthHitData> readTruth(const std::string& inputDir, size_t event) {
        auto truth = readEverything<ActsExamples::TruthHitData>(inputDir, "truth.csv", {}, event);
        return truth;
    }

    boost::unordered_map<uint64_t,std::vector<uint64_t>> readTruth_pid(const std::string& inputDir, size_t event) {
        auto truth = readEverythingStruct(inputDir, "truth.csv", {}, event);
        return truth;
    }

    std::vector<ActsExamples::ParticleData> readParticle(const std::string& inputDir, size_t event) {
        auto particle = readEverything<ActsExamples::ParticleData>(inputDir, "particles_initial.csv", {}, event);
        return particle;
    }

}  // namespace

ActsExamples::ProcessCode ActsExamples::CsvHitReader::read(
    const ActsExamples::AlgorithmContext& ctx) {
  // hit_id in the files is not required to be neither continuous nor
  // monotonic. internally, we want continous indices within [0,#hits)
  // to simplify data handling. to be able to perform this mapping we first
  // read all data into memory before converting to the internal event data
  // types.

  auto hitsData = readHits(m_cfg.inputDir, ctx.eventNumber);
  auto truthData = readTruth(m_cfg.inputDir, ctx.eventNumber);
  auto particleData = readParticle(m_cfg.inputDir, ctx.eventNumber);
  auto truthHit_pid = readTruth_pid(m_cfg.inputDir, ctx.eventNumber);

  // Prepare containers for the hit data using the framework event data types
  //typedef boost::bimap<Acts::GeometryIdentifier,struct Hitinformation> bimap;
  GeometryIdMultimap<struct Hitinformation> hits;
  typedef std::multimap<ActsFatras::Barcode,struct Hitinformation> map_hitsPID;
  map_hitsPID hitsPID;
  
  hits.reserve(hitsData.size());

  Hitinformation hitInfo;
  for (const HitData& hit : hitsData) {
    Acts::GeometryIdentifier::Value geoId = hit.geometry_id;
    hitInfo.geometry_id = geoId;
    hitInfo.x=hit.x * Acts::UnitConstants::mm;
    hitInfo.y=hit.y * Acts::UnitConstants::mm;
    hitInfo.z=hit.z * Acts::UnitConstants::mm;
    hitInfo.hit_id=hit.hit_id;


    //Get the particle_id and the index of the hit along the track
    auto truthRange = makeRange(std::equal_range(truthData.begin(), truthData.end(), hit.hit_id, CompareHitId{}));
    for (const auto& truth : truthRange) {
        hitInfo.particle_id=truth.particle_id;
        hitInfo.index=truth.index;
    }
    //following not well coded
    std::vector<Acts::GeometryIdentifier::Value> seen;
    for (auto& i : truthData){
        
        if ((hitInfo.particle_id == i.particle_id)){
            seen.push_back(i.geometry_id);
            if ((i.geometry_id == geoId) && (hitInfo.hit_id != i.hit_id)){
                hitInfo.duplicate = true;
                if (i.index> hitInfo.index) {
                    //std::cout<<"\n"<<hitInfo.hit_id<<std::endl;
                    hitInfo.linkUp=true;
                    //std::cout<<hitInfo.particle_id<<" "<<i.particle_id<<std::endl;
                    //std::cout<<geoId<<" "<<i.geometry_id<<std::endl;
                    //std::cout<<hitInfo.index<<" "<<i.index<<std::endl;
                    //std::cout<<hitInfo.hit_id<<" "<<i.hit_id<<std::endl;
                    //std::cout<<hitInfo.duplicate<<" "<<hitInfo.linkUp<<std::endl;
                    break;
                }else{
                    hitInfo.linkUp=false;
                }
            }else{
            hitInfo.duplicate = false; 
            }
        }
    }
    hitInfo.nhits_wihoutDuplicate=0;
    
    for(size_t i=0; i<seen.size();i++){
        for (size_t y=0; y<seen.size();y++){
            if (seen[i]==seen[y] && i!=y){
                hitInfo.nhits_wihoutDuplicate=-1;
            }
        }
    }
    

    ActsFatras::Barcode barcode = ActsFatras::Barcode(hitInfo.particle_id);
    //get the pt of the particle (initial)
    auto particleRange = makeRange(std::equal_range(particleData.begin(), particleData.end(),
                                        hitInfo.particle_id, CompareParticleId{}));
    for (const auto& p : particleRange) hitInfo.pt_particle=sqrt(p.px*p.px + p.py*p.py) * Acts::UnitConstants::GeV;  
    

    //Get the number of hits of this particle (useful to remove unwanted true connection)
    auto particleTrueRange = makeRange(truthHit_pid.equal_range(hitInfo.particle_id));
    for(auto& p: particleTrueRange) hitInfo.nhit_particle = p.second.size();
    
    hitInfo.nhits_wihoutDuplicate = hitInfo.nhits_wihoutDuplicate + hitInfo.nhit_particle;

    hits.insert(std::pair<Acts::GeometryIdentifier,ActsExamples::CsvHitReader::Hitinformation>(geoId,hitInfo));
    hitsPID.insert(std::pair<ActsFatras::Barcode,ActsExamples::CsvHitReader::Hitinformation>(barcode,hitInfo));
  }

  // Write the data to the EventStore
  ctx.eventStore.add(m_cfg.outputHits_map, std::move(hits));
  ctx.eventStore.add(m_cfg.outputHits_vec, std::move(hitsPID));

  return ActsExamples::ProcessCode::SUCCESS;
}
