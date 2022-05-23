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


#include "thz-spectrum-propagation-loss.h"
#include <ns3/thz-spectrum-waveform.h>
#include <ns3/mobility-model.h>
#include <ns3/object.h>
#include <ns3/double.h>
#include <ns3/antenna-model.h>
#include <ns3/cosine-antenna-model.h>
#include <ns3/angles.h>
#include <ns3/core-module.h>
#include <ns3/log.h>
#include <cmath> // for M_PI
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <random>
#include <algorithm>
#include <unistd.h>
#include <stdio.h>
#include <limits.h>

NS_LOG_COMPONENT_DEFINE ("THzSpectrumPropagationLoss");

namespace ns3 {

NS_OBJECT_ENSURE_REGISTERED (THzSpectrumPropagationLoss);

THzSpectrumPropagationLoss::THzSpectrumPropagationLoss ()
{
}

THzSpectrumPropagationLoss::~THzSpectrumPropagationLoss ()
{
}

       
Ptr<SpectrumValue>
THzSpectrumPropagationLoss::CalcRxPowerSpectralDensity (Ptr<const SpectrumValue> txPsd,
                                                        Ptr<const MobilityModel> a,
                                                        Ptr<const MobilityModel> b) 
{
  Ptr<SpectrumValue> rxPsd = Copy<SpectrumValue> (txPsd); //[W]
  Values::iterator vit = rxPsd->ValuesBegin ();
  Bands::const_iterator fit = rxPsd->ConstBandsBegin ();

  NS_ASSERT (a);
  NS_ASSERT (b);
  double d = a->GetDistanceFrom (b);
  while (vit != rxPsd->ValuesEnd ())
    {
      NS_ASSERT (fit != rxPsd->ConstBandsEnd ());
      double lossDb = 10 * std::log10 (CalculateSpreadLoss (fit->fc, d)) + 10 * std::log10 (CalculateAbsLoss (fit->fc, d));
      double lossW = std::pow (10.0, lossDb / 10.0);
      *vit /= lossW;
      ++vit;
      ++fit;
    }
  return rxPsd; //[W]
}

double
THzSpectrumPropagationLoss::CalcRxPowerDA (Ptr<THzSpectrumSignalParameters> txParams,
                                           Ptr<MobilityModel> a,
                                           Ptr<MobilityModel> b,
                                           double RxTxGainDb) 
{
  double RxTxGainW = std::pow (10.0, (RxTxGainDb) / 10.0);
  Ptr<SpectrumValue> rxPsd = Copy<SpectrumValue> (txParams->txPsd); // [W]
  Values::iterator vit = rxPsd->ValuesBegin ();
  Bands::const_iterator fit = rxPsd->ConstBandsBegin ();

  NS_ASSERT (a);
  NS_ASSERT (b);
  double d = a->GetDistanceFrom (b),i = 0;
  NS_LOG_INFO ("Distance = " << d);
  double rxPsd_inte = 0.0;
  while (vit != rxPsd->ValuesEnd ())
    {
      NS_ASSERT (fit != rxPsd->ConstBandsEnd ());
      double lossDb = 10 * std::log10 (CalculateSpreadLoss (fit->fc, d)) + 10 * std::log10 (CalculateAbsLoss (fit->fc, d));
      double lossW = std::pow (10.0, lossDb / 10.0);
      *vit /= lossW;
      rxPsd_inte += *vit;
      ++vit;
      ++fit;
      ++i;
    }
  NS_LOG_INFO ("rxPsd_inte = " << rxPsd_inte << " RxTxGainW " << RxTxGainW);
  double rxPower = rxPsd_inte * txParams->subBandBandwidth * (txParams->numberOfSubBands / txParams->numberOfSamples) * RxTxGainW;
  double rxPowerDbm = 10 * std::log10 (rxPower * 1000.0);
  NS_LOG_INFO ("Number of samples: " << txParams->numberOfSamples << " RxPower = " << rxPower << " W ");
  NS_LOG_INFO ("RxPowerDbm = " << rxPowerDbm << " Dbm ");
  return rxPowerDbm;
}


//-----------------------------------------------------------------------------------------------------//
double
THzSpectrumPropagationLoss::CalculateSpreadLoss (double f, double d) const
{
  NS_ASSERT (d >= 0);

  if (d == 0)
    {
      return 0;
    }

  NS_ASSERT (f > 0);
  double loss_sqrt = (4 * M_PI * f * d) / 3e8;
  double loss = loss_sqrt * loss_sqrt;
  return loss;
}

double
THzSpectrumPropagationLoss::CalculateAbsLoss (double f, double d) 
{
  double kf = 0.0;
  double loss = 0.0;

  if(mapContainsKey(m_freqMap, f))
  {
    kf = m_freqMap[f];
  }
  else
  {
    std::ifstream AbsCoefile;
    AbsCoefile.open ("contrib/thz/model/data_AbsCoe.txt", std::ifstream::in);
    if (!AbsCoefile.is_open ())
        {
        NS_FATAL_ERROR ("THzSpectrumPropagationLoss::CalculateAbsLoss: open data_AbsCoe.txt failed 1");
        }

    std:: ifstream frequencyfile;
    frequencyfile.open ("contrib/thz/model/data_frequency.txt", std::ifstream::in);
    if (!frequencyfile.is_open ())
        {
        NS_FATAL_ERROR ("THzSpectrumPropagationLoss::CalculateAbsLoss: open data_frequency.txt failed");
        }
    double f_ite;
    double k_ite;
    int i = 0;
    int j = 0;

    while (frequencyfile >> f_ite)
        {
        if (f_ite < f - 9.894e8 || f_ite > f + 9.894e8)
            {
            j++;
            }
        else
            {
            break;
            }
        }
    while (AbsCoefile >> k_ite)
        {
        if (i != j)
            {
            i++;
            }
        else
            {
            kf = k_ite;
            break;
            }
        NS_ASSERT (d >= 0);

        if (d == 0)
            {
            return 0;
            }
        }
    NS_ASSERT (f > 0);
    m_freqMap.insert(std::pair<double, double> (f, kf));
    NS_LOG_UNCOND("inserted to map f: " << f << " kf: "<< kf);
  } // if f != m_previousFc

  loss = exp (kf * d);
  return loss;
}


Ptr<SpectrumValue>
THzSpectrumPropagationLoss::LoadedAbsCoe (int s, int j, double f, double d,Ptr<const SpectrumValue> txPsd) const
{
  std:: ifstream AbsCoefile;
  AbsCoefile.open ("contrib/thz/model/data_AbsCoe.txt", std::ifstream::in);
  if (!AbsCoefile.is_open ())
    {
      NS_FATAL_ERROR ("THzSpectrumPropagationLoss::LoadedAbsCoe: open data_AbsCoe.txt failed");
    }
  double k;
  Ptr<SpectrumValue> kf_store = Copy <SpectrumValue> (txPsd);
  int i = 1;

  while (AbsCoefile >> k)
    {
      if (i < s)
        {
          i++;
        }
      else if (i > j)
        {
          break;
        }
      else
        {
          (*kf_store)[i - s] = k;
          i++;
        }
    }
  return kf_store;
}

bool 
THzSpectrumPropagationLoss::mapContainsKey(std::map<double, double>& map, double key)
{
  if (map.find(key) == map.end()) return false;
  return true;
}

double
THzSpectrumPropagationLoss::DoCalcHybridModelRxPower ( Ptr<const THzSpectrumSignalParameters> txParams,
                                                       Ptr<const MobilityModel> aMob,
                                                       Ptr<const MobilityModel> bMob,
                                                       Ptr<const THzDirectionalAntenna> aAntenna,
                                                       Ptr<const THzDirectionalAntenna> bAntenna,
                                                       const uint16_t maxRay,
                                                       const double rxTxGainDb,
                                                       const std::string m_qdPath) const
{
  NS_LOG_FUNCTION (this);
  
  Ptr<const MatrixBasedChannelModel::ChannelMatrix> channelParamH = GetHBChannel (aMob, bMob, aAntenna, bAntenna, m_qdPath);
  double rxPowerDbm = RxPowerCal (channelParamH, maxRay, rxTxGainDb, txParams);

  return rxPowerDbm;
}


double
THzSpectrumPropagationLoss::DoCalcFSModelRxPower (Ptr<const THzSpectrumSignalParameters> txParams,
                                                  Ptr<const MobilityModel> aMob,
                                                  Ptr<const MobilityModel> bMob,
                                                  Ptr<const THzDirectionalAntenna> aAntenna,
                                                  Ptr<const THzDirectionalAntenna> bAntenna,
                                                  const uint16_t maxRay,
                                                  const double rxTxGainDb) const
{
  NS_LOG_FUNCTION (this);
  
  Ptr<const MatrixBasedChannelModel::ChannelMatrix> channelParamS = GetFullyStatChannel (aMob, bMob, aAntenna, bAntenna);
  double rxPowerDbm = RxPowerCal (channelParamS, maxRay, rxTxGainDb, txParams);

  return rxPowerDbm;
}


Ptr<const MatrixBasedChannelModel::ChannelMatrix>
THzSpectrumPropagationLoss::GetHBChannel (Ptr<const MobilityModel> aMob,
                                          Ptr<const MobilityModel> bMob,
                                          Ptr<const THzDirectionalAntenna> aAntenna,
                                          Ptr<const THzDirectionalAntenna> bAntenna,
                                          const std::string m_qdPath) const
{
  NS_LOG_FUNCTION (this);

  // instance qd channel and stat channels
  // HB = qd + stat
  Ptr< MatrixBasedChannelModel::ChannelMatrix> qdChannel = GetRTChannel (aMob, bMob, aAntenna, bAntenna, m_qdPath);
  Ptr< const MatrixBasedChannelModel::ChannelMatrix> statChannel = GetStatChannel (qdChannel, aMob, bMob, aAntenna, bAntenna);
  
  qdChannel->m_delay.insert (qdChannel->m_delay.end (), statChannel->m_delay.begin (), statChannel->m_delay.end ());
  qdChannel->m_angle.insert (qdChannel->m_angle.end (), statChannel->m_angle.begin (), statChannel->m_angle.end ());
  for (uint8_t i = 0; i < statChannel->m_channel[0][0].size (); i++)
  {
    qdChannel->m_channel[0][0].push_back (statChannel->m_channel[0][0][i]);
  }

  for (uint8_t i = 0; i < qdChannel->m_delay.size (); i++)
  {
    NS_LOG_DEBUG("CHannel created with this params:" <<
                "Delay:" << qdChannel->m_delay[i]<<
                ",angle [0]:" << qdChannel->m_angle[0][i] <<
                ",angle [1]:" << qdChannel  ->m_angle[1][i] <<
                ",angle [2]:" << qdChannel->m_angle[2][i] <<
                ",angle [3]:" << qdChannel->m_angle[3][i] <<
                "channel:" <<qdChannel->m_channel[0][0][i]);
  }

  return qdChannel;
}


Ptr< MatrixBasedChannelModel::ChannelMatrix>
THzSpectrumPropagationLoss::GetRTChannel (Ptr<const MobilityModel> aMob,
                                          Ptr<const MobilityModel> bMob,
                                          Ptr<const THzDirectionalAntenna> aAntenna,
                                          Ptr<const THzDirectionalAntenna> bAntenna,
                                          const std::string m_qdPath)const
{
  NS_LOG_FUNCTION (this);

  // channelParams indicates prime channel model
  // channelParamsN indicates qd channel included antenna pattern
  Ptr<const MatrixBasedChannelModel::ChannelMatrix> channelParams = Create< MatrixBasedChannelModel::ChannelMatrix> ();
  Ptr< MatrixBasedChannelModel::ChannelMatrix> channelParamsN = Create< MatrixBasedChannelModel::ChannelMatrix> ();
  
  std::string qdFilesPath = "contrib/qd-channel/model/QD/"; // The path of the folder with the QD scenarios
 
  static Ptr<QdChannelModel> qdChannel = CreateObject<QdChannelModel> (qdFilesPath, m_qdPath);
  
  channelParams = qdChannel->GetRTChannel (aMob, bMob, aAntenna);

  channelParamsN->m_delay.clear ();
  channelParamsN->m_angle.clear ();
  channelParamsN->m_channel.clear ();
  std::copy (channelParams->m_delay.begin (), channelParams->m_delay.end (), std::back_inserter (channelParamsN->m_delay));
  std::copy (channelParams->m_angle.begin (), channelParams->m_angle.end (), std::back_inserter (channelParamsN->m_angle));
  std::copy (channelParams->m_channel.begin (), channelParams->m_channel.end (), std::back_inserter (channelParamsN->m_channel));
  channelParamsN->m_generatedTime = channelParams->m_generatedTime;
  channelParamsN->m_nodeIds = channelParams->m_nodeIds;

  double mpcCount = channelParams->m_channel[0][0].size ();
  NS_LOG_DEBUG ("mpcCount:" << mpcCount);

 std::complex<double>  txAntennaPattern;
 std::complex<double>  rxAntennaPattern;
 for (double mpcIndex = 1; mpcIndex < mpcCount; mpcIndex++)  // pass the '0' Mpc --> los scenarios antenna pattern included in the power calc step
 {
    txAntennaPattern = std::pow (10, (aAntenna->GetAntennaGainDbAngle (DegreesToRadians(channelParamsN->m_angle[2][mpcIndex]), aMob, bMob) / 20.0));
    rxAntennaPattern = std::pow (10, (aAntenna->GetAntennaGainDbAngle (DegreesToRadians(channelParamsN->m_angle[0][mpcIndex]), aMob, bMob) / 20.0));
    NS_LOG_DEBUG ("Antenna Gains are TX:" << txAntennaPattern << ", RX:" << rxAntennaPattern << ", total antenna gain is:" << txAntennaPattern * rxAntennaPattern);
    double delay = -2 * M_PI * 140e9 * (channelParamsN->m_delay[mpcIndex]);
    channelParamsN->m_channel[0][0][mpcIndex] = channelParamsN->m_channel[0][0][mpcIndex] * txAntennaPattern * rxAntennaPattern * std::exp (std::complex<double> (0, delay))* std::exp (std::complex<double> (0, channelParamsN->m_angle[0][mpcIndex]));
 }
    
  return channelParamsN;
}


Ptr<const MatrixBasedChannelModel::ChannelMatrix>
THzSpectrumPropagationLoss::GetStatChannel (Ptr<const MatrixBasedChannelModel::ChannelMatrix> qdChannel,
                                            Ptr<const MobilityModel> aMob,
                                            Ptr<const MobilityModel> bMob,
                                            Ptr<const THzDirectionalAntenna> aAntenna,
                                            Ptr<const THzDirectionalAntenna> bAntenna)const
{
  NS_LOG_FUNCTION (this);

  Ptr<MatrixBasedChannelModel::ChannelMatrix> channelParams = Create<MatrixBasedChannelModel::ChannelMatrix> ();

  NS_ASSERT (aMob);
  NS_ASSERT (bMob);
  uint32_t aId = aMob->GetObject<Node> ()->GetId ();
  uint32_t bId = bMob->GetObject<Node> ()->GetId ();

  double d = bMob->GetDistanceFrom (aMob);
  NS_LOG_DEBUG("aID:" << aId <<
               ",bID:"<< bId <<
               ",Distance:" << d);

  // Lrt is number RT clusters and ls is non ray tracing clusters using eq.14
  uint16_t LrtCount = qdChannel->m_channel[0][0].size () - 1; // los has no cluster
  uint16_t LSCount = ceil (-1.056 * d + 9.776);

  // SubPath Counts for pre and post of LT and LS Rays
  uint16_t RTPreSC = (int) LogNormalDist (2.09, 1.20);
  uint16_t RTPostSC = (int) LogNormalDist (4, 1.59);
  uint16_t LSPreSC = (int) LogNormalDist (1.76, 1.33);
  uint16_t LSPostSC = (int) LogNormalDist (2.07, 1.03);
  uint16_t RTTotSC = RTPreSC + RTPostSC;
  uint16_t LSTotSC = LSPreSC + LSPostSC + 1; // central path is included
  
  // total rays in Statisctical Channel
  uint16_t TotalRayCount = LrtCount * RTTotSC + LSCount * LSTotSC;
  
  NS_LOG_DEBUG(" RT count:"<< LrtCount <<
               ",LS count:"<< LSCount <<
               ",RTPreSC:"<< RTPreSC <<
               ",RTPostSC:"<< RTPostSC <<
               ",LSPreSC:"<< LSPreSC <<
               ",LSPostSC:"<< LSPostSC <<
               ",RTTotSC:"<< RTTotSC <<
               ",LSTotSC:"<< LSTotSC <<
               ",TotalRayCount:"<< TotalRayCount);
  

  // Channel Matrix params
  MatrixBasedChannelModel::Complex3DVector H;
  MatrixBasedChannelModel::Double2DVector sAngle;
  std::vector<double> sDelay;

  sDelay.resize (TotalRayCount);


  H.resize (1);
  for (uint64_t bIndex = 0; bIndex < 1; bIndex++)
  {
    H[bIndex].resize (1);
    for (uint64_t aIndex = 0; aIndex < 1; aIndex++)
    {
      if (TotalRayCount > 0)
      {
        H[bIndex][aIndex].resize (TotalRayCount, std::complex<double> (0, 0));
      }
      else
      {
        NS_LOG_DEBUG("No ray excists and no Matrix Genearated!!");
        H[bIndex][aIndex].resize (0, std::complex<double> (0, 0));
      }
    }
  }
	
  sAngle.resize (4);
  for (uint64_t aIndex = 0; aIndex < 4; aIndex++)
  {
    if (TotalRayCount > 0)
    {
      sAngle[aIndex].resize (TotalRayCount, 0);
    }
    else
    {
      NS_LOG_DEBUG ("No ray excists and no Matrix Genearated!!");
      sAngle[aIndex].resize (0, 0);
    }
  }

  // Param Inititlize, Values are updates due to TABLE VI of "channel measurment and Ray-tracing-statistical Hybrid modeling for Low-terahert Indoor Communication"
  double mu = 0;
  double kappa = 0;
  double landa = 0;
  double a = 0;
  double b = 0;
  std::complex<double> pathLoss;
  std::complex<double> ray;
  double plLos = std::norm (qdChannel->m_channel[0][0][0]);
  double mpcIndex = 0;
  double delay = 0;


	// Channel structure
	//  '1'                               '2'           '3'
	// RT {intra Cluster Rays} + RS { main Ray + intra Cluster Rays}	
	

	// '1' Generate RT  intra cluster Rays
	for (uint8_t i = 0; i < LrtCount; i++)
  {
    for (uint8_t j = 0; j < RTTotSC; j++)
    {   

      if (j < RTPreSC)
      {
        landa = 0.0918;
        mu = 0;
        kappa = 33;
        a = 0.41;
        b = 0.51;
      }
      else
      {
        landa = 0.0593;
        mu = 0.5;
        kappa = 3.5; 
        a = 0.42;
        b = 0.45; 
      }

      // Channel Calc For each MPC
      if (j == 0) // Calc Delay for the first ray of each cluster from the central ray of same cluster in RT part
      { 
        sDelay[mpcIndex] = qdChannel->m_delay[i + 1] + 1e-9 * ExpDist (landa);
      } 
      else
      {
        sDelay[mpcIndex] = sDelay[mpcIndex - 1] + 1e-9 * ExpDist (landa);        
      }

      sAngle[0][mpcIndex] = VonMisesDist (mu, kappa);		//AOA 
      sAngle[1][mpcIndex] = VonMisesDist (mu, kappa);		//ZOA 
      sAngle[2][mpcIndex] = qdChannel->m_angle[2][i+1];
      sAngle[3][mpcIndex] = qdChannel->m_angle[3][i+1];
      delay = -2 * M_PI * 140e9 * (sDelay[mpcIndex]);
      ray = std::norm (qdChannel->m_channel[0][0][i + 1]) * std::exp (std::complex<double> (0, delay)) * std::exp (std::complex<double> (0, sAngle[0][mpcIndex]));
      H[0][0][mpcIndex] += ray;
      mpcIndex++;
    }
  }	
		
  
  // '2' central ray of non-RT(LS) cluster geneartion

  // Params for central rays of clusters
	mu = -1.23;
	kappa = 0; 
	landa = 13.2;
	a = 0.36;
	b = 0.25;

  for (uint8_t i = 0; i < LSCount; i++)
  {
    if (i == 0) //First cluster in LS depends on Last Cluster of RT
    { 
      sDelay[mpcIndex] = sDelay[mpcIndex-1] + 1e-9 * ExpDist (landa);
    }
    else
    {
	    sDelay[mpcIndex] = sDelay[mpcIndex-LSTotSC]+ 1e-9 * ExpDist (landa);
	  }

	  sAngle[0][mpcIndex] = VonMisesDist(mu, kappa);		//AOA
    sAngle[1][mpcIndex] = VonMisesDist(mu, kappa);		//ZOA
    if (sAngle[1][mpcIndex] < 0)
      {
       sAngle[1][mpcIndex] += M_PI; 
      } 
	  pathLoss = SubPathLossCalc (plLos, a, b, qdChannel->m_delay[0], sDelay[mpcIndex]);
    delay = -2 * M_PI * 140e9 * (sDelay[mpcIndex]);
    ray =  pathLoss * std::exp (std::complex<double> (0, delay)) * std::exp (std::complex<double> (0, sAngle[0][mpcIndex]));
    H[0][0][mpcIndex] += ray;
    mpcIndex += LSTotSC;
  }


  mpcIndex = LrtCount * RTTotSC;  
	// '3' Generate other rays of non-RT(LS) cluster Rays
	for (uint8_t i = 0; i < LSCount; i++)
  {
    mpcIndex++; //  pass first central ray of clusters
    for (uint8_t j = 1; j< LSTotSC; j++)
    {
      if (j<LSPreSC)
      {
        landa = 0.102;
		    mu = 0;
        kappa = 33;
		    a = 0.41;
        b = 0.51;
      }
      else
      {
        landa =  0.1417;
		    mu = 0.5;
        kappa = 3.5; 
		    a = 0.42;
        b = 0.45; 
      }
      // Channel Calc For each MPC	
      sDelay[mpcIndex] = sDelay[mpcIndex - 1] + 1e-9 * ExpDist (landa);
      sAngle[0][mpcIndex] = VonMisesDist (mu, kappa);		//AOA 
      sAngle[1][mpcIndex] = VonMisesDist (mu, kappa);		//ZOA

      pathLoss = SubPathLossCalc (plLos, a, b, qdChannel->m_delay[0], sDelay[mpcIndex]);
      delay = -2 * M_PI * 140e9 * (sDelay[mpcIndex]);
      ray =  pathLoss * std::exp (std::complex<double> (0, delay)) * std::exp (std::complex<double> (0, sAngle[0][mpcIndex]));
      H[0][0][mpcIndex] += ray;
      mpcIndex++;
    }
  }


  channelParams->m_channel = H;

  channelParams->m_angle.clear ();
  channelParams->m_angle = sAngle;

  channelParams->m_delay.clear ();
  channelParams->m_delay = sDelay;

  channelParams->m_generatedTime = Simulator::Now ();

  channelParams->m_nodeIds = std::make_pair (aId, bId);
  NS_LOG_INFO ("STAT Channel generated, total rays in Stat Channel is:" << channelParams->m_channel[0][0].size ());

  return channelParams;
}

Ptr<const MatrixBasedChannelModel::ChannelMatrix>
THzSpectrumPropagationLoss::GetFullyStatChannel (Ptr<const MobilityModel> aMob,
                                                 Ptr<const MobilityModel> bMob,
                                                 Ptr<const THzDirectionalAntenna> aAntenna,
                                                 Ptr<const THzDirectionalAntenna> bAntenna) const
{
  NS_LOG_FUNCTION(this);
  Ptr<MatrixBasedChannelModel::ChannelMatrix> channelParams = Create<MatrixBasedChannelModel::ChannelMatrix> ();

  
  NS_ASSERT (aMob);
  NS_ASSERT (bMob);
  uint32_t aId = aMob->GetObject<Node> ()->GetId ();
  uint32_t bId = bMob->GetObject<Node> ()->GetId ();

  double d = bMob->GetDistanceFrom (aMob);
  // NS_LOG_DEBUG("aID:" << aId <<
  //              ",bID:"<< bId <<
  //              ",Distance:" << d);

  // Number of Time Cluster for LOS & NLOS
  uint16_t tcLos = PoissionDist (0.9) + 1;
  uint16_t tcNLos = PoissionDist (1.8) + 1;
 // double TC = tcLos +tcNLos;    //Total Time Cluster


  // Number of SubPaths in the clusters
  uint16_t tcLosSub = int (ExpDist (1.4)) + 1;
  uint16_t tcNLosSub = int (ExpDist (1.2)) + 1;
  NS_LOG_DEBUG("num of tcLos:" << tcLosSub <<
              ",num of tcNLos:" << tcNLosSub);
 // double TCSub = tcLosSub + tcNLosSub;

//Scenario LOS or NLOS
 Ptr<UniformRandomVariable> x = CreateObject<UniformRandomVariable> ();
  double value = x->GetValue (0, 1);
  bool scenario = false;
  if (value >= 0.7)
  {
    scenario = true; //LOS
  }

 uint16_t TotalRayCount = 0;
  if (scenario == true)
  {
    TotalRayCount = tcLos * tcLosSub;
  }
  else 
  {
    TotalRayCount = tcNLos * tcNLosSub; 
  }
  // Number of total Rays
  
  NS_ABORT_MSG_IF (TotalRayCount == 0, "NO Ray founded!");

  // NS_LOG_DEBUG("tcLos:"<< tcLos <<
  //             ",tcNLos:"<< tcNLos <<
  //             ",TC:" << TC << 
  //             ",tcLossSub:" << tcLosSub <<
  //             ",tcNLosSub:" << tcNLosSub <<
  //             ",TCSub:" << TCSub <<
  //             ", Total Rays Counts:" << TotalRayCount);

  // Channel Matrix Vectors
  MatrixBasedChannelModel::Complex3DVector H;
  MatrixBasedChannelModel::Double2DVector fsAngle;
  std::vector<double> fsDelay;
  std::vector<double> fsPhase;

  fsDelay.resize(TotalRayCount);
  fsPhase.resize(TotalRayCount);
  H.resize (1);

  for (uint64_t bIndex = 0; bIndex < 1; bIndex++)
  {
    H[bIndex].resize (1);
    for (uint64_t aIndex = 0; aIndex < 1; aIndex++)
    {
      if (TotalRayCount > 0)
      {
        H[bIndex][aIndex].resize (TotalRayCount, std::complex<double> (0, 0));
      } 
      else
      {
        H[bIndex][aIndex].resize (0, std::complex<double> (0, 0));
      }
    }
  }
	
  fsAngle.resize (4);
  for (uint64_t aIndex = 0; aIndex < 4; aIndex++)
  {
    if (TotalRayCount > 0)
    {
      fsAngle[aIndex].resize (TotalRayCount, 0);
    }
    else
    {
      NS_LOG_DEBUG("No ray excists and no Matrix Generated!!");
      fsAngle[aIndex].resize (0, 0);
      }
  }

	// Param Inititlize
  double MTI = 6e-9;    // 6 ns time Interval in Indoor Office - 140GHz refer to FS chnnel model paper
  double landa = 0;   // Delay interval
  double txAntennaPattern = 1;
  double rxAntennaPattern = 1;
  std::complex<double> ray = 0;
  double mpcIndex = 0;
  double lAod = 0;
  double lAoa = 0;
  double plCluster = 0;
  double plSub = 0;
  double delay = 0;

  // LOS Params
  double lAodMax = 2;
  double lAoaMax = 2;
  double thetaAod = 0;
  double thetaAoa = 0;
  double phiZod = 0;
  double phiZoa = 0;
  double deltaZod = 3.4;
  double deltaZoa = 3.3;
  double deltaAoa = 4.4;
  double deltaAod = 4.3;
  double ple = 2.1;
  double sigma = 9.1;
  double sigmaS = 4.6;
  double delayConst = 18.2;
  double delayConstS = 2;


if (scenario == true)
{
	// LOS MPCs
	for (uint8_t i = 0; i < tcLos; i++) //Clusters
  {
    lAod = DUniformDist (1, lAodMax);
    lAoa = DUniformDist (1, lAoaMax);
    thetaAod = UniformDist (360 * (lAod - 1) / lAodMax, 360 * (lAod / lAodMax));
    thetaAoa = UniformDist (360 * (lAoa - 1) / lAodMax, 360 * (lAoa / lAoaMax));        
    phiZod = NormDist (-6.8, 4.9);
    phiZoa = NormDist (7.4, 4.5);

    for (uint8_t j = 0; j < tcLosSub; j++)  // Sub Rays in Clusters
    {
      if (j == 0) // Central Ray
      {
        landa = 14.6;
        if (i == 0) //first ray Delay 
        {
          fsDelay[0] = MTI + ExpDist (landa) * 1e-9; 
        }
        else
        {
          fsDelay[mpcIndex] = fsDelay[mpcIndex - 1] + MTI + 1e-9 * ExpDist (landa);
        }
        // Path loss of each cluster calcs onces in each cluster
        plCluster = FSClusterPathLossCalcdB (d, ple, sigma, fsDelay[mpcIndex], delayConst);
      }
      else 
      {
        landa = 1.1;
        fsDelay[mpcIndex] = fsDelay[mpcIndex - 1] + 1e-9 * ExpDist (landa);
      }
      
      fsPhase[mpcIndex] = UniformDist (0, 2 * M_PI);   //RAD
      fsAngle[0][mpcIndex] = NormDist (thetaAoa, deltaAoa); 
      fsAngle[1][mpcIndex] = NormDist (phiZoa, deltaZoa);
      fsAngle[2][mpcIndex] = NormDist (thetaAod, deltaAod); 
      fsAngle[3][mpcIndex] = NormDist (phiZod, deltaZod);

      // NS_LOG_DEBUG("LOS,,,lAod:"<< lAod <<
      //               ",lAoa:"<< lAoa <<
      //               ",thetaAod:" << thetaAod << 
      //               ",thetaAoa:" << thetaAoa <<
      //               ",fsAngle[0][mpcIndex]:" << fsAngle[0][mpcIndex] <<
      //               ",fsAngle[2][mpcIndex]:" << fsAngle[2][mpcIndex]);

      plSub = std::pow (10, (FSSubPathLossCalcdB (plCluster, sigmaS, fsDelay[mpcIndex], delayConstS)) / 20); 
      txAntennaPattern = pow (10,(aAntenna->GetAntennaGainDbAngle (DegreesToRadians(fsAngle[2][mpcIndex]), aMob, bMob) / 20.0));
      rxAntennaPattern = pow (10,(bAntenna->GetAntennaGainDbAngle (DegreesToRadians(fsAngle[0][mpcIndex]), aMob, bMob) / 20.0));
      delay = -2 * M_PI * 140e9 * (fsDelay[mpcIndex]);

      if (i == 0 && j == 0) // antenna Gain applied in power calc step for LOS
      {
        ray = plSub * exp (std::complex<double> (0, delay));
      }
      else
      {
        ray = plSub * txAntennaPattern * rxAntennaPattern * exp (std::complex<double> (0, delay));
      }
      H[0][0][mpcIndex] += ray;
      mpcIndex++;
    }
  }	
}
else 
{	
	// NLOS Params
  lAodMax = 2;
  lAoaMax = 2;
  thetaAod = 0;
  thetaAoa = 0;
  phiZod = 0;
  phiZoa = 0;
  deltaZod = 3.3;
  deltaZoa = 3.3;
  deltaAoa = 5.6;
  deltaAod = 4.0;
  ple = 4.6;
  sigma = 12.8;
  sigmaS = 5.8;
  delayConst = 16.1;
  delayConstS = 2.4;

  // Generate Rays for NLOS
	for (uint8_t i = 0; i < tcNLos; i++)
  {
    lAod = DUniformDist (1, lAodMax);
    lAoa = DUniformDist (1, lAoaMax);
    thetaAod = UniformDist (360 * (lAod - 1) / lAodMax, 360 * (lAod/ lAodMax));
    thetaAoa = UniformDist (360 * (lAoa - 1) / lAodMax, 360 * (lAoa/ lAoaMax));
    phiZod = NormDist (-2.5, 2.7);
    phiZoa = NormDist (4.8, 2.8);

    for (uint8_t j = 0; j < tcNLosSub; j++)
    {
      if (j == 0) // Central Ray
      {
        landa = 21;
        fsDelay[mpcIndex] = fsDelay[mpcIndex - 1] + MTI + 1e-9 * ExpDist (landa);
        // Path loss of each cluster calcs onces in each cluster
        plCluster = FSClusterPathLossCalcdB (d, ple, sigma, fsDelay[mpcIndex], delayConst);
      }
      else 
      {
        landa = 2.7;
        fsDelay[mpcIndex] = fsDelay[mpcIndex - 1] + 1e-9 * ExpDist (landa);
      }
      fsPhase[mpcIndex] = UniformDist (0, 2 * M_PI);   //RAD
      fsAngle[0][mpcIndex] = NormDist (thetaAoa, deltaAoa); 
      fsAngle[1][mpcIndex] = NormDist (phiZoa, deltaZoa);
      fsAngle[2][mpcIndex] = NormDist (thetaAod, deltaAod); 
      fsAngle[3][mpcIndex] = NormDist (phiZod, deltaZod);
        // NS_LOG_DEBUG("NLOS,,,lAod:"<< lAod <<
        //               ",lAoa:"<< lAoa <<
        //               ",thetaAod:" << thetaAod << 
        //               ",thetaAoa:" << thetaAoa <<
        //               ",fsAngle[0][mpcIndex]:" << fsAngle[0][mpcIndex] <<
        //               ",fsAngle[2][mpcIndex]:" << fsAngle[2][mpcIndex]);

      plSub = std::pow (10, (FSSubPathLossCalcdB (plCluster, sigmaS, fsDelay[mpcIndex], delayConstS)) / 20);
      txAntennaPattern = std::pow (10, (aAntenna->GetAntennaGainDbAngle (DegreesToRadians (fsAngle[2][mpcIndex]), aMob, bMob) / 20.0));
      rxAntennaPattern = std::pow (10, (bAntenna->GetAntennaGainDbAngle (DegreesToRadians (fsAngle[0][mpcIndex]), aMob, bMob) / 20.0));
      delay = -2 * M_PI * 140e9 * (fsDelay[mpcIndex]);
      ray = plSub * txAntennaPattern * rxAntennaPattern * exp (std::complex<double> (0, delay));
      H[0][0][mpcIndex] += ray;
      mpcIndex++;
    }
  }	
}  
  channelParams->m_channel = H;

  channelParams->m_angle.clear();
  channelParams->m_angle = fsAngle;

  channelParams->m_delay.clear();
  channelParams->m_delay = fsDelay;

  channelParams->m_generatedTime = Simulator::Now ();

  channelParams->m_nodeIds = std::make_pair (aId, bId);

  NS_LOG_INFO("Fully STAT Channel generated, total rays is:" <<channelParams->m_channel[0][0].size());

  for (uint8_t i = 0; i< channelParams->m_delay.size(); i++)
  {
    NS_LOG_DEBUG("CHannel created with this params:" <<
                ", Delay:" << channelParams->m_delay[i]<<
                ",angle [0]:" << channelParams->m_angle[0][i] <<
                ",angle [1]:" << channelParams->m_angle[1][i] <<
                ",angle [2]:" << channelParams->m_angle[2][i] <<
                ",angle [3]:" << channelParams->m_angle[3][i] <<
                "channel:" << channelParams->m_channel[0][0][i]);

  }

  return channelParams;
}                                               


double
THzSpectrumPropagationLoss::FSClusterPathLossCalcdB(const double distance, const double ple, const double sigma, const double delay, const double coefDelay) const
{

  double c = 299792458;
  double fc = 142e9;
  double FSPL = 20 * log10 (4 * M_PI * fc / c);
  double shadow = LogNormalDist (0, sigma);

  double PL = FSPL + 10 * ple * log10 (distance) + shadow;

  double clusterPLdB = PL + 10 * log10 (exp (-delay / coefDelay));

  NS_LOG_DEBUG("Params in Path Loss Calc, distance:" << distance <<
              ",shadow:" << shadow <<
              ",FSPL:" <<FSPL<<
              ",ple:" << ple<<
              ",PL:" <<PL <<
              ",Cluster PL:" <<clusterPLdB);

  return clusterPLdB;
}


double
THzSpectrumPropagationLoss::FSSubPathLossCalcdB (double clusterPL, const double sigma, const double delay, const double coefDelay) const
{
  double shadow = LogNormalDist (0, sigma);
  double clusterPLdB = clusterPL + 10 * log10 (exp (-delay / coefDelay)) + shadow;

 NS_LOG_DEBUG("Value of power (dB) is FS Channel is:" << clusterPLdB <<
              ", main val is (dB):" << clusterPL);
  return -1 * clusterPLdB;
}



double 
THzSpectrumPropagationLoss::PoissionDist (double landaC) const
{

  NS_LOG_FUNCTION (this << landaC);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::poisson_distribution<int> distribution(landaC);
  double value = distribution (gen);
 // NS_LOG_DEBUG(" Valus in Pois Dist is:" << value);
  return value;
 }


double 
THzSpectrumPropagationLoss::NormDist (const double mean, const double variance) const
{
  //NS_LOG_FUNCTION (this << mu << sigma);
  Ptr<NormalRandomVariable> x = CreateObject<NormalRandomVariable> ();

  x->SetAttribute ("Mean", DoubleValue (mean));
  x->SetAttribute ("Variance", DoubleValue (variance));
 // x->SetAttribute ("Bound", DoubleValue (bound));

  double value = x->GetValue ();

  return value;
 }


double 
THzSpectrumPropagationLoss::UniformDist (const double min, const double max) const
{
  NS_LOG_FUNCTION (this << min << max);
  Ptr<UniformRandomVariable> x = CreateObject<UniformRandomVariable> ();
  double value = x->GetValue (min, max);
  NS_LOG_DEBUG(" Valus in Uniform Dist is:" << (double) value <<",min is:"<<min<<", max is:"<<max);
  return value;
 }

uint16_t 
THzSpectrumPropagationLoss::DUniformDist (const double min, const double max) const
{

  NS_LOG_FUNCTION (this << min << max);
  Ptr<UniformRandomVariable> x = CreateObject<UniformRandomVariable> ();
  uint16_t value = x->GetInteger (min, max);
 // NS_LOG_DEBUG(" Valus in Uniform Dist is:" << (double) value <<",min is:"<<min<<", max is:"<<max);
  return value;
 }


double 
THzSpectrumPropagationLoss::ExpDist (const double landa) const
{
  NS_LOG_FUNCTION (this);
  Ptr<ExponentialRandomVariable> x = CreateObject<ExponentialRandomVariable> ();
  x->SetAttribute ("Mean", DoubleValue (landa));
  double value = x->GetValue();
  //NS_LOG_DEBUG("ExpDist :"<< landa << ",:" <<value);
  return value;
}

double 
THzSpectrumPropagationLoss::LogNormalDist (const double mu, const double sigma) const
{
  //NS_LOG_FUNCTION (this << mu << sigma);

  Ptr<LogNormalRandomVariable> x = CreateObject<LogNormalRandomVariable> ();

  x->SetAttribute ("Mu", DoubleValue (mu));
  x->SetAttribute ("Sigma", DoubleValue (sigma));

  double value = x->GetValue ();
  double quntized = std::min(value, mu + sigma);
  NS_LOG_DEBUG("ns log dist :"<< value << ", quntize:" <<quntized);

  return quntized;
 }

double 
THzSpectrumPropagationLoss::VonMisesDist (const double mu, const double kappa) const
{
//  NS_LOG_FUNCTION (this << mu << kappa);
  double theta = 0;
  if (kappa == 0)
  {
  Ptr<UniformRandomVariable> x = CreateObject<UniformRandomVariable> ();
  theta= x->GetValue (-M_PI, M_PI);
  }

  else
  {
  Ptr<NormalRandomVariable> x = CreateObject<NormalRandomVariable> ();
  x->SetAttribute ("Mean", DoubleValue (0));
  x->SetAttribute ("Variance", DoubleValue (1/kappa));
  theta =  x->GetValue() * M_PI;  // Random value -PI::PI
  }

  return theta;
}



std::complex<double> 
THzSpectrumPropagationLoss::SubPathLossCalc (const double pathLoss, const double a, const double b, const double TOALos, const double TOAPath) const
{
  NS_LOG_FUNCTION (this << pathLoss << a << b << TOALos << TOAPath);
  //NS_LOG_DEBUG("VAL inside of Func is:"<<pathLoss);
  std::complex<double> lossVal = 0;
  lossVal = pathLoss * a *(pow(abs((TOALos-TOAPath) * 1e9), b));
  NS_LOG_DEBUG("Path Loss VAL inside of Func is:"<<lossVal <<", PathLosss was:"<<pathLoss);
  return lossVal;
}


double 
THzSpectrumPropagationLoss::RxPowerCal (const Ptr<const MatrixBasedChannelModel::ChannelMatrix> channelparamS,
                                        const uint16_t maxRay,
                                        const double rxTxGainDb,
                                        const Ptr<const THzSpectrumSignalParameters> txParams) const
{
  uint16_t numCluster = channelparamS->m_channel[0][0].size ();
  double activeRays = std::min (numCluster, maxRay);
  NS_LOG_INFO("numCluster:" << numCluster);
  double RxTxGainW = std::pow (10.0, (rxTxGainDb) / 20.0);
  
  Ptr<SpectrumValue> rxPsd = Copy<SpectrumValue> (txParams->txPsd); // [W]
  Values::iterator vit = rxPsd->ValuesBegin ();
  Bands::const_iterator fit = rxPsd->ConstBandsBegin ();
  double rxPsd_inte = 0.0;
  double gainW = 0;

  while (vit != rxPsd->ValuesEnd ())
    {
      NS_LOG_DEBUG("txPsd is:" << *vit);
      NS_ASSERT (fit != rxPsd->ConstBandsEnd ());

      gainW = (std::norm (channelparamS->m_channel[0][0][0]) * RxTxGainW);  //LOS line
      for (uint8_t cIndex = 1; cIndex < activeRays; cIndex++)
      {
        gainW += (std::norm(channelparamS->m_channel[0][0][cIndex]));
        NS_LOG_UNCOND("value of gain W is: " << gainW << " at time step of: " << Now());
      }
        NS_LOG_DEBUG("gain is (w):" <<gainW);
      *vit *= gainW;
      rxPsd_inte += *vit;
      ++vit;
      ++fit;
    }
  NS_LOG_INFO ("rxPsd_inte = " << rxPsd_inte);
  double rxPower = rxPsd_inte * txParams->subBandBandwidth * (txParams->numberOfSubBands / txParams->numberOfSamples);
  double rxPowerDbm = 10 * std::log10 (rxPower * 1000.0);
  NS_LOG_INFO ("Number of samples: " << txParams->numberOfSamples << " RxPower = " << rxPower << " W ");
  NS_LOG_INFO ("RxPowerDbm = " << rxPowerDbm << " Dbm ");
  
  return rxPowerDbm;
}


} // namespace ns3
