// -----------------------------------------------------------------------------
//
// Copyright (C) The BioDynaMo Project.
// All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------

#ifndef CORE_EXECUTION_CONTEXT_IN_PLACE_EXEC_CTXT_H_
#define CORE_EXECUTION_CONTEXT_IN_PLACE_EXEC_CTXT_H_

#include <functional>
#include <memory>
#include <utility>
#include <vector>

#include "core/container/so_uid_map.h"
#include "core/functor.h"
#include "core/operation/operation.h"
#include "core/sim_object/so_uid.h"
#include "core/util/spinlock.h"
#include "core/util/thread_info.h"

namespace bdm {

class SimObject;

/// This execution context updates simulation objects in place. \n
/// Let's assume we have two sim objects `A, B` in our simulation that we want
/// to update to the next timestep `A*, B*`. If we have one thread it will first
/// update `A` and afterwards `B` and write the updates directly to the same
/// data structure. Therefore, before we start updating `B` the array looks
/// like this: `A*, B`. `B` already observes the updated `A`. \n
/// Operations in method `Execute` are executed in order given by the user.
/// Subsequent operations observe the changes of earlier operations.\n
/// In-place updates can lead to race conditions if simulation objects not only
/// modify themselves, but also neighbors. Therefore, a protection mechanism has
/// been added. \see `Param::thread_safety_mechanism_`
/// New sim objects will only be visible at the next iteration. \n
/// Also removal of a sim object happens at the end of each iteration.
class InPlaceExecutionContext {
 public:
  struct ThreadSafeSoUidMap {
    using value_type = std::pair<SimObject*, uint64_t>;
    ThreadSafeSoUidMap();
    ~ThreadSafeSoUidMap();

    void Insert(const SoUid& uid, const value_type& value);
    const value_type& operator[](const SoUid& key);
    uint64_t Size() const;
    void Resize(uint64_t new_size);
    void RemoveOldCopies();

    using Map = SoUidMap<value_type>;
    Spinlock lock_;
    Spinlock next_lock_;
    Map* map_;
    Map* next_;
    std::vector<Map*> previous_maps_;
  };

  explicit InPlaceExecutionContext(
      const std::shared_ptr<ThreadSafeSoUidMap>& map);

  virtual ~InPlaceExecutionContext();

  /// This function is called at the beginning of each iteration to setup all
  /// execution contexts.
  /// This function is not thread-safe.
  /// NB: Invalidates references and pointers to simulation objects.
  void SetupIterationAll(
      const std::vector<InPlaceExecutionContext*>& all_exec_ctxts) const;

  /// This function is called at the end of each iteration to tear down all
  /// execution contexts.
  /// This function is not thread-safe. \n
  /// NB: Invalidates references and pointers to simulation objects.
  void TearDownIterationAll(
      const std::vector<InPlaceExecutionContext*>& all_exec_ctxts) const;

  /// Execute a series of operations on a simulation object in the order given
  /// in the argument
  void Execute(SimObject* so, const std::vector<Operation*>& operations);

  void push_back(SimObject* new_so);  // NOLINT

  void ForEachNeighbor(Functor<void, const SimObject*, double>& lambda,
                       const SimObject& query);

  /// Forwards the call to `Grid::ForEachNeighborWithinRadius`
  void ForEachNeighborWithinRadius(
      Functor<void, const SimObject*, double>& lambda, const SimObject& query,
      double squared_radius);

  SimObject* GetSimObject(const SoUid& uid);

  const SimObject* GetConstSimObject(const SoUid& uid);

  void RemoveFromSimulation(const SoUid& uid);

 private:
  /// Lookup table SoUid -> SoPointer for new created sim objects
  std::shared_ptr<ThreadSafeSoUidMap> new_so_map_;

  ThreadInfo* tinfo_;

  /// Contains unique ids of sim objects that will be removed at the end of each
  /// iteration.
  std::vector<SoUid> remove_;
  std::vector<Spinlock*> locks;

  /// Pointer to new sim objects
  std::vector<SimObject*> new_sim_objects_;

  /// prevent race conditions for cached SimObjects
  std::atomic_flag mutex_ = ATOMIC_FLAG_INIT;

  std::vector<std::pair<const SimObject*, double>> neighbor_cache_;
};

}  // namespace bdm

#endif  // CORE_EXECUTION_CONTEXT_IN_PLACE_EXEC_CTXT_H_
