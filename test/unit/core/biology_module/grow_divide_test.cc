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

#include "core/biology_module/grow_divide.h"
#include <typeinfo>
#include "core/resource_manager.h"
#include "gtest/gtest.h"
#include "unit/test_util/test_util.h"

namespace bdm {

TEST(GrowDivideTest, Grow) {
  Simulation simulation(TEST_NAME);
  auto* rm = simulation.GetResourceManager();
  auto* ctxt = simulation.GetExecutionContext();

  ctxt->SetupIterationAll(simulation.GetAllExecCtxts());

  Cell cell;
  cell.SetDiameter(40);

  GrowDivide gd(40, 300, {gAllEventIds});
  gd.Run(&cell);

  ctxt->TearDownIterationAll(simulation.GetAllExecCtxts());

  EXPECT_NEAR(33513.321638291127, cell.GetVolume(), abs_error<double>::value);
  EXPECT_EQ(0u, rm->GetNumSimObjects());
}

TEST(GrowDivideTest, Divide) {
  Simulation simulation(TEST_NAME);
  auto* rm = simulation.GetResourceManager();
  auto* ctxt = simulation.GetExecutionContext();

  ctxt->SetupIterationAll(simulation.GetAllExecCtxts());

  Cell cell;
  cell.SetDiameter(41);

  GrowDivide gd(40, 300, {gAllEventIds});
  gd.Run(&cell);

  ctxt->TearDownIterationAll(simulation.GetAllExecCtxts());

  EXPECT_GT(41, cell.GetDiameter());
  EXPECT_EQ(1u, rm->GetNumSimObjects());
}

}  // namespace bdm
