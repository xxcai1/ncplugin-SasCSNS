#include "NCSansIsotropic.hh"
//Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCRandUtils.hh"
#include <vector>



NCP::SansIsotropic::SansIsotropic( const std::vector<double>& Q, const std::vector<double>& intensity )
: m_iofq(Q, intensity), 
  m_xs(std::unique_ptr<LookUpTable>(nullptr))
{
  //Parse and validate values:
  if(Q.size()!=intensity.size())
    NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                     <<" (Q.size()!=intensity.size()" );  
  // total xs
  std::vector<double> envec, xsvec;
  envec.reserve(Q.size());
  xsvec.reserve(Q.size());


  // neutron integral cross section is
  // sigma(E) = 2*pi/k_o^2 \int_0^{2*k_o} Q*I(Q) dQ
  // For a pointwise data, we can use the trapezoidal rule to integrate.
  // We first creat an energy grid E correspond to the Q grid, each E_i corresponds to the phase space 
  // that covers by the I(Q) in the range of between Q_0 and Q_i, .i.e E[0]=NC::k2ekin(0.5*Q[0]). So,
  // sigma(E_{i+1}) = 2*pi /k_o^2 \sum 0.5 * [Q_{i+1} - I(Q_i)] [Q_{i+1} * I(Q_{i+1}) + Q_i * I(Q_i)]
  //                = pi /k_o^2 \sum [Q_{i+1} - I(Q_i)] [Q_{i+1} * I(Q_{i+1}) + Q_i * I(Q_i)]
  // when i=0, Sigma(E1) = 4*pi*I(Q).
  

  double en=NC::wl2ekin(4*NC::kPi/Q[0]);
  envec.push_back(en);
  xsvec.push_back(4*NC::kPi*intensity[0]);
  double accumIntegrand(0.);

  for(size_t i=1;i<Q.size();i++)
  {
    double lastK=2.*NC::kPi/NC::ekin2wl(envec.back()); //ekin to ki
    accumIntegrand = xsvec.back()*lastK*lastK;
    accumIntegrand += NC::kPi*(Q[i]-Q[i-1])*(intensity[i]*Q[i]+intensity[i-1]*Q[i-1]);
    double k = Q[i]*0.5;
    envec.push_back(NC::wl2ekin(2*NC::kPi/k));
    xsvec.push_back(accumIntegrand/(k*k));
  }
  m_xs=std::make_unique<LookUpTable>(envec, xsvec, NCP::LookUpTable::Extrapolate::kConst_OverSqrtX);

}

double NCP::SansIsotropic::calcCrossSection( double neutron_ekin ) const
{
  if(m_xs.get())
    return m_xs->get(neutron_ekin);
  else
    return 0.;
}

NCP::SansIsotropic::ScatEvent NCP::SansIsotropic::sampleScatteringEvent( NC::RNG& rng, double neutron_ekin ) const
{
  ScatEvent result;
  result.ekin_final = neutron_ekin;
  double Q = m_iofq.sampleQValue(rng, NC::NeutronEnergy(neutron_ekin));
  double kappa = NC::ekin2k(neutron_ekin);
  result.mu = 1-Q*Q/(2*kappa*kappa);
  return result;
}
