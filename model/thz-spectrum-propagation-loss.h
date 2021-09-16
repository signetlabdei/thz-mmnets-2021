/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2021 Northeastern University (https://unlab.tech/)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * Author: Qing Xia <qingxia@buffalo.edu>
 *         Zahed Hossain <zahedhos@buffalo.edu>
 *         Josep Miquel Jornet <j.jornet@northeastern.edu>
 *         Daniel Morales <danimoralesbrotons@gmail.com>
 */


#ifndef THZ_SPECTRUM_PROPAGATION_LOSS_H
#define THZ_SPECTRUM_PROPAGATION_LOSS_H

#include <ns3/object.h>
#include <map>
#include <ns3/mobility-model.h>
#include <ns3/spectrum-value.h>
#include "thz-spectrum-signal-parameters.h"
#include "ns3/matrix-based-channel-model.h"
#include "ns3/thz-dir-antenna.h"
#include "ns3/qd-channel-utils.h"
#include "ns3/qd-channel-model.h"

namespace ns3 {
/**
 * \defgroup Terahertz Spectrum Propagation Loss Models
 *
 */

class THzSpectrumPropagationLoss : public Object
{
public:
  THzSpectrumPropagationLoss ();
  virtual ~THzSpectrumPropagationLoss ();

  virtual bool mapContainsKey(std::map<double, double>& map, double key);

  /**
    * \param txPsd the power spectral density of the transmitted signal, unit in Watt.
    * \param a the mobility of sender.
    * \param b the mobility of receiver.
    *
    * \return the power spectral density of the received signal, unit is Watt.
    *
    * This application doesn't apply the terahertz directional antenna.
    */
  virtual Ptr<SpectrumValue> CalcRxPowerSpectralDensity (Ptr<const SpectrumValue> txPsd,
                                                         Ptr<const MobilityModel> a,
                                                         Ptr<const MobilityModel> b) ;

  /**
    * \param txPsd the power spectral density of the transmitted signal, unit in Watt.
    * \param a the mobility of sender.
    * \param b the mobility of receiver.
    * \param RxTxGainDb the total antenna gain of both transmitter and receiver, unit in dB.
    *
    * \return the received signal power, unit in dBm.
    *
    * This application applies the terahertz directional antenna.
    */
  virtual double CalcRxPowerDA (Ptr<THzSpectrumSignalParameters> txParams,
                                Ptr<MobilityModel> a,
                                Ptr<MobilityModel> b,
                                double RxTxGainDb) ;

  /**
    * \brief Calculate the spreading loss
    *
    * \param f the central frequency of the operation subband, unit in Hz.
    * \param d the distance between transmitter and receiver, unit in meter.
    *
    * \return the spreading loss unit in Watt.
    *
    * The reference is J. M. Jornet and I. F. Akyildiz, Channel modeling and capacity analysis
    * of electromagnetic wireless nanonetworks in the terahertz band, IEEE
    * Transactions on Wireless Communications, vol. 10, no. 10, pp. 3211
    * 3221, Oct. 2011
    *
    * The values of f and d are collected from HITRAN database.
    */
  virtual double CalculateSpreadLoss (double f, double d) const;

  /**
    * \brief Calculate the absorption coefficient loss.
    *
    * \param f the central frequency of the operation subband, unit in Hz.
    * \param d the distance between transmitter and receiver, unit in meter.
    *
    * \return the absorption coefficient loss, unit in Watt.
    *
    * The reference is J. M. Jornet and I. F. Akyildiz, Channel modeling and capacity analysis
    * of electromagnetic wireless nanonetworks in the terahertz band, IEEE
    * Transactions on Wireless Communications, vol. 10, no. 10, pp. 3211
    * 3221, Oct. 2011.
    *
    * The values of f and d are collected from HITRAN database.
    */
  virtual double CalculateAbsLoss (double f, double d) ;

  /**
    * \param s the starting boundary of the absorption coefficent.
    * \param j the ending boundary of the absorption coefficent.
    * \param f the central frequency of the operation subband, unit in Hz.
    * \param d the distance between transmitter and receiver, unit in meter.
    * \param txPsd the power spectral density of the transmitted signal.
    *
    * \return the list of absorption coefficent values within the specified boundaries.
    *
    * Mainly for checking purpose.
    */
  virtual Ptr<SpectrumValue> LoadedAbsCoe (int s, int j, double f, double d,Ptr<const SpectrumValue> txPsd) const;

  double m_previousFc;
  double m_kf;
  std::map<double, double> m_freqMap;

   /**
    * \param txPsd the power spectral density of the transmitted signal, unit in Watt.
    * \param aMob the mobility of sender.
    * \param bMob the mobility of receiver.
    * \param aAntenna the antenna model of sender.
    * \param bAntenna the antenna model of receiver.
    * \param maxRay the maximum number of rays that can be counted to calculate the power.
    * \param RxTxGainDb the total antenna gain of both transmitter and receiver, unit in dB.
    * \param m_qdPath path to the qd channel model params from Ray tracing 
    *
    * \return the power spectral density of the received signal, unit is Watt.
    *
    */
  double DoCalcHybridModelRxPower (Ptr<const THzSpectrumSignalParameters> txParams,
                                   Ptr<const MobilityModel> aMob,
                                   Ptr<const MobilityModel> bMob,
                                   Ptr<const THzDirectionalAntenna> aAntenna,
                                   Ptr<const THzDirectionalAntenna> bAntenna,
                                   const uint16_t maxRay,
                                   const double rxTxGainDb,
                                   const std::string m_qdPath ) const;
  /**
    * \param txPsd the power spectral density of the transmitted signal, unit in Watt.
    * \param aMob the mobility of sender.
    * \param bMob the mobility of receiver.
    * \param aAntenna the antenna model of sender.
    * \param bAntenna the antenna model of receiver.
    * \param maxRay the maximum number of rays that can be counted to calculate the power.
    * \param RxTxGainDb the total antenna gain of both transmitter and receiver, unit in dB.
    *
    * \return the power spectral density of the received signal, unit is Watt.
    *
    */
  double DoCalcFSModelRxPower (Ptr<const THzSpectrumSignalParameters> txParams,
                               Ptr<const MobilityModel> aMob,
                               Ptr<const MobilityModel> bMob,
                               Ptr<const THzDirectionalAntenna> aAntenna,
                               Ptr<const THzDirectionalAntenna> bAntenna,
                               const uint16_t maxRay,
                               const double rxTxGainDb) const;

  /**
    * \param aMob the mobility of sender.
    * \param bMob the mobility of receiver.
    * \param aAntenna the antenna model of sender.
    * \param bAntenna the antenna model of receiver.
    * \param m_qdPath path to the qd channel model params from Ray tracing 
    * \return the Hybrid channel model include rays, pathloss, angles, and delays.
    *
    */
  Ptr<const MatrixBasedChannelModel::ChannelMatrix> GetHBChannel (Ptr<const MobilityModel> aMob,
                                                                  Ptr<const MobilityModel> bMob,
                                                                  Ptr<const THzDirectionalAntenna> aAntenna,
                                                                  Ptr<const THzDirectionalAntenna> bAntenna,
                                                                  const std::string m_qdPath) const;

  /**
    * \param aMob the mobility of sender.
    * \param bMob the mobility of receiver.
    * \param aAntenna the antenna model of sender.
    * \param bAntenna the antenna model of receiver.
    * \param m_qdPath path to the qd channel model params from Ray tracing 
    *
    * \return the channel model of Ray tracing include rays, pathloss, angles, and delays.
    *
    */
  Ptr< const MatrixBasedChannelModel::ChannelMatrix> GetRTChannel (Ptr<const MobilityModel> aMob,
                                                                   Ptr<const MobilityModel> bMob,
                                                                   Ptr<const THzDirectionalAntenna> aAntenna,
                                                                   Ptr<const THzDirectionalAntenna> bAntenna,
                                                                   const std::string m_qdPath) const;

  /**
    * \param aMob the mobility of sender.
    * \param bMob the mobility of receiver.
    * \param aAntenna the antenna model of sender.
    * \param bAntenna the antenna model of receiver.
    *
    * \return the statistical part of Hybrid channel model include rays, pathloss, angles, and delays.
    *
    */
  Ptr<const MatrixBasedChannelModel::ChannelMatrix> GetStatChannel (Ptr<const MatrixBasedChannelModel::ChannelMatrix> qdchannel,
                                                                    Ptr<const MobilityModel> aMob,
                                                                    Ptr<const MobilityModel> bMob,
                                                                    Ptr<const THzDirectionalAntenna> aAntenna,
                                                                    Ptr<const THzDirectionalAntenna> bAntenna) const;

  /**
    * \param aMob the mobility of sender.
    * \param bMob the mobility of receiver.
    * \param aAntenna the antenna model of sender.
    * \param bAntenna the antenna model of receiver.
    *
    * \return the fully statisctical channel model include rays, pathloss, angles, and delays.
    *
    */
  Ptr<const MatrixBasedChannelModel::ChannelMatrix> GetFullyStatChannel (Ptr<const MobilityModel> aMob,
                                                                         Ptr<const MobilityModel> bMob,
                                                                         Ptr<const THzDirectionalAntenna> aAntenna,
                                                                         Ptr<const THzDirectionalAntenna> bAntenna) const;

   /**
    * \param distance  distance between transmitter and receiver, unit in meter.
    * \param ple the value of ple for LOS and NLOS.
    * \param sigma the sigma for distribution of power.
    * \param delay delay of the ray in nanoseconds.
    * \param coefDelay the coefficiente for calculating the delay effect in power.
    *
    * \return the power in dB.
    *
    */
  double FSClusterPathLossCalcdB (const double distance, const double ple, const double sigma, const double delay, const double coefDelay)const;

   /**
    * \param clusterPL  power of the cluster in db.
    * \param sigma the sigma for distribution of power.
    * \param delay delay of the ray in nanoseconds.
    * \param coefDelay the coefficiente for calculating the delay effect in power.
    *
    * \return the power in dB.
    *
    */
  double FSSubPathLossCalcdB(double clusterPL, const double sigma, const double delay, const double coefDelay)const;

   /**
    * \param landaC the coefficent for calculating poission distribution.
    *
    * \return the probability based on the possion dist and landac
    *
    */
  double PoissionDist(const double landaC) const;

    /**
    * \param mean the coefficent declares the mean of distribution.
    * \param var the coefficent declares the variance of distribution.
    *
    * \return the probability
    *
    */
  double NormDist(const double mean, const double var) const;

    /**
    * \param min the coefficent declares the min value of distribution.
    * \param max the coefficent declares the max value of distribution.
    *
    * \return the probability in double
    *
    */
  double UniformDist(const double min, const double max) const;

    /**
    * \param min the coefficent declares the min value of distribution.
    * \param max the coefficent declares the max value of distribution.
    *
    * \return the probability in int
    *
    */
  uint16_t DUniformDist(const double min, const double max) const;

    /**
    * \param mu the coefficent in exponential distribution.
    *
    * \return the probability in double
    *
    */
  double ExpDist (const double mu) const;

    /**
    * \param mean the coefficent declares the mean of distribution.
    * \param bound the coefficent declares the maximum bound of return value in distribution.
    *
    * \return the probability in double
    *
    */
  double LogNormalDist (const double mean, const double bound) const;

    /**
    * \param mu the coefficent declares the measure of the location of distribution.
    * \param kappa the coefficent declares  the measure of the concentration of distribution.
    *
    * \return the probability in double
    *
    */
  double VonMisesDist (const double mu, const double kappa) const;

   /**
    * \param pathLoss the starting boundary of the absorption coefficent.
    * \param a the coefficent for calculating power.
    * \param b other coefficent for calculating power related to Hbrid channel model.
    * \param ToaLos the delay of ray to reach the reciver.
    * \param ToaPath delay of the main path.
    *
    * \return the power value in dB.
    *
    * Mainly for checking purpose.
    */
  std::complex<double> SubPathLossCalc (const double pathLoss, const double a, const double b, const double ToaLos, const double ToaPath) const;


   /**
    * \param channelparamS the channel model include rays, pathloss, angles, and delays.
    * \param maxRay the maximum number of rays that can be counted to calculate the power.
    * \param rxTxGainDb the total antenna gain of both transmitter and receiver, unit in dB.
    * \param txPsd the power spectral density of the transmitted signal, unit in Watt.
    *
    * \return the power value in dB.
    *
    * Mainly for checking purpose.
    */
  double RxPowerCal(const Ptr<const MatrixBasedChannelModel::ChannelMatrix> channelparamS,
                    const uint16_t maxRay,
                    const double rxTxGainDb,
                    const Ptr<const THzSpectrumSignalParameters> txParams) const;
  
};


} // namespace ns3

#endif /* SPECTRUM_PROPAGATION_LOSS_H */


