#ifndef PARAM_H_
#define PARAM_H_

namespace bdm {

class Param {
 public:
  /// Time between two simulation step, in hours.
  static constexpr float kSimulationTimeStep = 0.01;
  /// Maximum jump that a point mass can do in one time step. Useful to
  /// stabilize the simulation
  static constexpr float kSimulationMaximalDisplacement = 3.0;

  static constexpr float kDefaultTolerance = 0.000000001;

  /// Maximum length of a discrete segment before it is cut into two parts.
  static constexpr float kNeuriteMaxLength = 20;  // usual value : 20
  /// Minimum length of a discrete segment before. If smaller it will try to
  /// fuse with the proximal one
  static constexpr float kNeuriteMinLength = 2.0;  // usual value : 10

  // Diffusion (saving time by not running diffusion for too small differences)

  /// If concentration of a substance smaller than this, it is not diffused
  static constexpr float kMinimalConcentrationForExtracellularDiffusion = 1e-5;

  /// If absolute concentration difference is smaller than that, there is no
  /// diffusion
  static constexpr float
      kMinimalDifferenceConcentrationForExtracacellularDiffusion = 1e-5;
  /// If ratio (absolute concentration difference)/concentration is smaller than
  /// that, no diffusion.
  static constexpr float kMinimalDCOverCForExtracellularDiffusion = 1e-3;

  /// If concentration of a substance smaller than this, it is not diffused
  static constexpr float kMinimalConcentrationForIntracellularDiffusion =
      1e-10;
  /// If absolute concentration difference is smaller than that, there is no
  /// diffusion
  static constexpr float
      kMinimalDifferenceConcentrationForIntracellularDiffusion = 1e-7;
  /// If ration (absolute concentration difference)/concentration is smaller
  /// than that, no diffusion.
  static constexpr float kMinimalDCOverCForIntracellularDiffusion = 1e-4;

  // Neurites
  /// Initial value of the restingLength before any specification.
  static constexpr float kNeuriteDefaultActualLength = 1.0;
  /// Diameter of an unspecified (= axon/dendrite) neurite when extends from the
  /// somaElement
  static constexpr float kNeuriteDefaultDiameter = 1.0;
  static constexpr float kNeuriteMinimalBifurcationLength = 0;
  /// Spring constant
  static constexpr float kNeuriteDefaultSpringConstant = 10;  // 10;
  /// Threshold the force acting on a neurite has to reach before a move is made
  /// ( = static friction).
  static constexpr float kNeuriteDefaultAdherence = 0.1;
  /// Rest to the movement ( = kinetic friction).
  static constexpr float kNeuriteDefaultMass = 1;

  static constexpr float kNeuriteDefaultTension = 0.0;

  // Somata
  /// CAUTION: not the radius but the diameter
  static constexpr float kSphereDefaultDiameter = 20;
  /// Threshold the force acting on a somaElement has to reach before a move is
  /// made ( = static friction).
  static constexpr float kSphereDefaultAdherence = 0.4;
  /// Restistance to the movement ( = kinetic friction).
  static constexpr float kSphereDefaultMass = 1;

  static constexpr float kSphereDefaultRotationalInertia = 0.5;

  static constexpr float kSphereDefaultInterObjectCoefficient = 0.15;

  /// Helpful constant to compare with 0
  static constexpr float kEpsilon = 1e-20;
  /// Helpful constant to identify 'infinity'
  static constexpr float kInfinity = 1e20;
};

}  // namespace bdm

#endif  // PARAM_H_
