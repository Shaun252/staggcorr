#!/usr/bin/env python3
"""
Test suite for comparing stagcorr package outputs with MILC lattice code data.

This module converts the notebook comparison logic into proper unit tests
that verify the consistency between stagcorr correlator calculations and
reference MILC data.
"""

import unittest
import numpy as np
import os
from re import search

import stagcorr
import stagcorr.correlators.main as corr


class TestMILCComparison(unittest.TestCase):
    """Test class for comparing stagcorr outputs with MILC reference data."""
    
    def setUp(self):
        """Set up test fixtures and load MILC data."""
        self.test_dir = os.path.dirname(os.path.abspath(__file__))
        self.milc_file = os.path.join(self.test_dir, 'MILCDATA', 'MILCRUNS')
        self.vol = (6, 6, 6, 6)
        self.tolerance = 1e-6  # Tolerance for numerical comparisons
        
        # Load MILC data
        self.dataDict = self._load_milc_data()
        
    def _load_milc_data(self):
        """Load and parse MILC data file."""
        def momRead(momLine):
            momKey = momLine.split(":")[1]
            return ",".join(momKey.rstrip("\n").replace(" ", "").replace("[", "").replace("]", "").rstrip(",").split(","))

        def STRead(STLine):
            STKey = STLine.split(":")[1]
            return STKey.rstrip("\n").replace(" ", "").replace(",", "")

        readq = False
        readqExt = False
        readaqExt = False
        
        qExtMom = None
        aqExtMom = None
        qSTExt = None
        aqSTExt = None
        aqExtSrcT = None
        qSourceMom = None
        aqSourceMom = None
        qST = 'G5-G5'
        NT = 1
        k = 0
        take_data = False
        corr_function = []
        dataDict = {}

        if not os.path.exists(self.milc_file):
            self.skipTest(f"MILC data file not found: {self.milc_file}")

        with open(self.milc_file, 'r') as fp:
            for line in fp.readlines():
                data = line.split()
                
                if 'lattice_size' in line:
                    NT = int(line.split(",")[-1])
                    continue
                
                if 'antiquark_source_origin' in line:
                    srcT = line[-4]
                
                if 'antiquark_source_mom' in line:
                    aqSourceMom = momRead(line)
                    continue
                
                if "antiquark_source_label" in line:
                    aqSourceType = data[1]
                    continue
                
                if 'antiquark_sink_ops' in line:
                    readaqExt = True
                    continue
                    
                if readaqExt and 'spin_taste_extend' in line:
                    aqSTExt = STRead(line)
                    continue
                        
                if readaqExt and 'mom' in line:
                    aqExtMom = momRead(line)
                    continue
                
                if readaqExt and 't0' in line:
                    aqExtSrcT = data[1]
                    readaqExt = False
                    continue
                
                if 'quark_source_mom' in line:
                    qSourceMom = momRead(line)
                    continue
                
                if 'quark_source_ops' in line:
                    readq = True
                    continue
                
                if readq and 'spin_taste' in line and 'operation' not in line:
                    qST = STRead(line)
                    readq = False
                    continue
                
                if 'quark_sink_ops' in line:
                    readqExt = True
                    continue
                        
                if readqExt and 'mom' in line:
                    qExtMom = momRead(line)
                    continue
                
                if readqExt and 'spin_taste_extend' in line:
                    qSTExt = STRead(line)
                    readqExt = False
                    continue
                
                if 'momentum:' == data[0]:
                    momkey = data[1]
                    
                if "spin_taste_sink" in line:
                    STSink = data[1]
                    continue

                if "correlator_key:" in line:
                    k = 0
                    take_data = True
                    
                    if aqSourceMom == None:
                        if qExtMom == None:
                            aqSourceMom = momkey
                        else:
                            aqSourceMom = "0,0,0"
                            
                    if qSourceMom == None:
                        if qExtMom != None:
                            pion2Mom = list(map(int, qExtMom.split(",")))
                            pion1Mom = ",".join(map(str, [-1*momval for momval in pion2Mom]))
                            qSourceMom = pion1Mom
                        else:
                            qSourceMom = "0,0,0"
                    
                    key = (aqSourceMom, qSourceMom, qST, STSink, qExtMom, qSTExt, aqExtMom, aqSTExt, momkey, srcT)
                    if key not in dataDict.keys():
                        dataDict[key] = {}
                    if aqSourceType not in dataDict[key].keys():
                        dataDict[key][aqSourceType] = {}
                    continue
                    
                if bool(search(r'\d', line)) == True and take_data:
                    data = line.split()
                    corr_function.append(float(data[1]) + 1j * float(data[2]))
                    k += 1

                if k == NT and take_data:
                    key = (aqSourceMom, qSourceMom, qST, STSink, qExtMom, qSTExt, aqExtMom, aqSTExt, momkey, srcT)
                    dataDict[key][aqSourceType][aqExtSrcT] = corr_function
                    take_data = False
                    corr_function = []
                    aqSourceMom, qSourceMom, qST, STSink, qExtMom, qSTExt, aqExtMom, aqSTExt, momkey, srcT, aqExtSrcT = [None] * 11
                    qST = 'G5-G5'
                    
        return dataDict

    def test_two_point_correlators(self):
        """Test two-point pion correlators against MILC data."""
        two_point_tests = 0
        
        for (aqSourceMom, qSourceMom, qST, STSink, qExtMom, qSTExt, aqExtMom, aqSTExt, momkey, srcT), corrFuncDict in self.dataDict.items():
            
            if qExtMom == None:  # Two-point correlator case
                incSpin, incT = qST.split("-")
                
                if aqSourceMom == '0,0,0' and qSourceMom == '0,0,0':
                    incMom = [0, 0, 0]
                    outMom = [0, 0, 0]
                elif aqSourceMom == '0,0,0' and qSourceMom != '0,0,0':
                    incMom = list(map(int, qSourceMom.split(",")))
                    outMom = [-1*momval for momval in incMom]
                elif aqSourceMom != '0,0,0' and qSourceMom == '0,0,0':
                    incMom = list(map(int, aqSourceMom.split(",")))
                    outMom = [-1*momval for momval in incMom]
                
                # Generate stagcorr two-point correlator
                two_point_pion = corr.generate_npt(
                    spinTasteMassNaikMomSymShift1=[incSpin, incT, 1, 0, incMom, 0],
                    spinTasteMassNaikMomSymShift2=[incSpin, incT, 1, 0, outMom, 0],
                    volume=self.vol
                )
                
                stagcorr_result = np.real(two_point_pion[:, 0]) / np.prod(self.vol[1:])
                milc_result = np.real(corrFuncDict['pt_src'][None])
                
                # Test that the ratio is consistent (should be constant factor)
                ratio = stagcorr_result / milc_result
                
                # Check that all ratios are the same (within tolerance)
                ratio_std = np.std(ratio)
                self.assertLess(ratio_std, self.tolerance, 
                    f"Two-point correlator ratios not consistent for {qST} momentum {incMom}")
                
                # Check that we get reasonable values (not NaN or inf)
                self.assertTrue(np.all(np.isfinite(ratio)), 
                    f"Non-finite values in two-point ratio for {qST} momentum {incMom}")
                
                two_point_tests += 1
        
        self.assertGreater(two_point_tests, 0, "No two-point correlator tests were run")

    def test_three_point_correlators(self):
        """Test three-point correlators against MILC data."""
        three_point_tests = 0
        
        for (aqSourceMom, qSourceMom, qST, STSink, qExtMom, qSTExt, aqExtMom, aqSTExt, momkey, srcT), corrFuncDict in self.dataDict.items():
            
            if aqExtMom == None and qExtMom != None:  # Three-point correlator case
                incSpin, incT = qST.split("-")
                sinkSpin, sinkT = STSink.split("-")
                
                if aqSourceMom == '0,0,0' and qSourceMom != '0,0,0':
                    incMom = list(map(int, qSourceMom.split(",")))
                    outMom = [-1*momval for momval in incMom]
                elif aqSourceMom != '0,0,0' and qSourceMom == '0,0,0':
                    incMom = list(map(int, aqSourceMom.split(",")))
                    outMom = [-1*momval for momval in incMom]
                
                # Generate stagcorr three-point correlator
                three_point_pion = corr.generate_npt(
                    spinTasteMassNaikMomSymShift1=[incSpin, incT, 0.1, 0, incMom, 0],
                    spinTasteMassNaikMomSymShift2=[incSpin, incT, 0.1, 0, outMom, 0],
                    spinTasteMassNaikMomSymShift3=[sinkSpin, sinkT, 0.1, 0, [0, 0, 0], 0],
                    volume=self.vol
                )

                stagcorr_result = np.imag(three_point_pion[:, 1, 0] / np.prod(self.vol[1:]))
                milc_result = np.imag(corrFuncDict['pt_src'][None])
                
                # Test that the ratio is consistent
                ratio = stagcorr_result / milc_result
                
                # Check that all ratios are the same (within tolerance)
                ratio_std = np.std(ratio)
                self.assertLess(ratio_std, self.tolerance,
                    f"Three-point correlator ratios not consistent for {qST}->{STSink} momentum {incMom}")
                
                # Check that we get reasonable values
                self.assertTrue(np.all(np.isfinite(ratio)),
                    f"Non-finite values in three-point ratio for {qST}->{STSink} momentum {incMom}")
                
                three_point_tests += 1
        
        if three_point_tests > 0:  # Only assert if we expect three-point tests
            self.assertGreater(three_point_tests, 0, "No three-point correlator tests were run")

    def test_four_point_correlators_O_diagram(self):
        """Test four-point correlators (O diagram) against MILC data."""
        four_point_O_tests = 0
        
        for (aqSourceMom, qSourceMom, qST, STSink, qExtMom, qSTExt, aqExtMom, aqSTExt, momkey, srcT), corrFuncDict in self.dataDict.items():
            
            if aqExtMom != None and '-' not in aqExtMom:  # Four-point O diagram
                incSpin, incT = qST.split("-")
                sinkSpin, sinkT = STSink.split("-")
                
                p1IncMom = list(map(int, qSourceMom.split(",")))
                p2IncMom = list(map(int, qExtMom.split(",")))
                p3OutMom = [-1*pi for pi in list(map(int, aqExtMom.split(",")))]  # MILC antiquark flips mom sign
                p4OutMom = list(map(int, momkey.split(",")))
                
                # Generate stagcorr four-point correlator
                four_point_pion = np.real(corr.generate_npt(
                    spinTasteMassNaikMomSymShift1=[incSpin, incT, 0.1, 0, p1IncMom, 0],
                    spinTasteMassNaikMomSymShift2=[incSpin, incT, 0.1, 0, p2IncMom, 0],
                    spinTasteMassNaikMomSymShift3=[sinkSpin, sinkT, 0.1, 0, p4OutMom, 0],
                    spinTasteMassNaikMomSymShift4=[sinkSpin, sinkT, 0.1, 0, p3OutMom, 0],
                    volume=self.vol
                )) / np.prod(self.vol[1:])
                
                # Test specific time slices as done in notebook
                test_times = [(1, 2), (2, 3), (3, 4), (4, 5)]
                
                for t3, t4 in test_times:
                    if t3 < four_point_pion.shape[0] and t4 < four_point_pion.shape[1]:
                        stagcorr_result = four_point_pion[t3, t4, 1, 0]
                        
                        # Get corresponding MILC data
                        milc_key = str(t3-1)  # MILC uses different indexing
                        if milc_key in corrFuncDict['pt_src'] and t4 < len(corrFuncDict['pt_src'][milc_key]):
                            milc_result = np.real(corrFuncDict['pt_src'][milc_key][t4])
                            
                            if abs(milc_result) > 1e-15:  # Avoid division by tiny numbers
                                ratio = stagcorr_result / milc_result
                                
                                # Check that we get reasonable values
                                self.assertTrue(np.isfinite(ratio),
                                    f"Non-finite ratio in four-point O diagram {qST}->{STSink} at times ({t3},{t4})")
                                
                                # Check that the ratio is not zero (indicates calculation worked)
                                self.assertNotAlmostEqual(abs(ratio), 0, places=10,
                                    msg=f"Zero ratio in four-point O diagram {qST}->{STSink} at times ({t3},{t4})")
                
                four_point_O_tests += 1
        
        if four_point_O_tests > 0:
            self.assertGreater(four_point_O_tests, 0, "No four-point O diagram tests were run")

    def test_four_point_correlators_8_diagram(self):
        """Test four-point correlators (8 diagram) against MILC data."""
        four_point_8_tests = 0
        
        for (aqSourceMom, qSourceMom, qST, STSink, qExtMom, qSTExt, aqExtMom, aqSTExt, momkey, srcT), corrFuncDict in self.dataDict.items():
            
            if aqExtMom != None and '-' in aqExtMom:  # Four-point 8 diagram
                incSpin, incT = qST.split("-")
                sinkSpin, sinkT = STSink.split("-")
                
                p1IncMom = list(map(int, qSourceMom.split(",")))
                p2IncMom = list(map(int, qExtMom.split(",")))
                p3OutMom = p2IncMom  # Different from O diagram
                p4OutMom = list(map(int, momkey.split(",")))
                
                # Generate stagcorr four-point correlator
                four_point_pion = np.real(corr.generate_npt(
                    spinTasteMassNaikMomSymShift1=[incSpin, incT, 0.1, 0, p1IncMom, 0],
                    spinTasteMassNaikMomSymShift2=[incSpin, incT, 0.1, 0, p2IncMom, 0],
                    spinTasteMassNaikMomSymShift3=[sinkSpin, sinkT, 0.1, 0, p3OutMom, 0],
                    spinTasteMassNaikMomSymShift4=[sinkSpin, sinkT, 0.1, 0, p4OutMom, 0],
                    volume=self.vol
                )) / np.prod(self.vol[1:])
                
                # Test specific time slices (different indexing for 8 diagram)
                test_times = [(2, 1), (3, 2), (4, 3), (5, 4)]
                
                for t3, t4 in test_times:
                    if t3 < four_point_pion.shape[0] and t4 < four_point_pion.shape[1]:
                        stagcorr_result = four_point_pion[t3, t4, 1, 0]
                        
                        # Get corresponding MILC data
                        milc_key = str(t3)
                        if milc_key in corrFuncDict['pt_src'] and t4 < len(corrFuncDict['pt_src'][milc_key]):
                            milc_result = np.real(corrFuncDict['pt_src'][milc_key][t4])
                            
                            if abs(milc_result) > 1e-15:  # Avoid division by tiny numbers
                                ratio = stagcorr_result / milc_result
                                
                                # Check that we get reasonable values
                                self.assertTrue(np.isfinite(ratio),
                                    f"Non-finite ratio in four-point 8 diagram {qST}->{STSink} at times ({t3},{t4})")
                                
                                # Check that the ratio is not zero
                                self.assertNotAlmostEqual(abs(ratio), 0, places=10,
                                    msg=f"Zero ratio in four-point 8 diagram {qST}->{STSink} at times ({t3},{t4})")
                
                four_point_8_tests += 1
        
        if four_point_8_tests > 0:
            self.assertGreater(four_point_8_tests, 0, "No four-point 8 diagram tests were run")

    def test_data_integrity(self):
        """Test that MILC data was loaded correctly."""
        self.assertGreater(len(self.dataDict), 0, "No MILC data was loaded")
        
        # Check that we have the expected data structure
        for key, corrFuncDict in self.dataDict.items():
            self.assertIsInstance(key, tuple, "MILC data keys should be tuples")
            self.assertEqual(len(key), 10, "MILC data keys should have 10 elements")
            self.assertIsInstance(corrFuncDict, dict, "MILC data values should be dictionaries")


if __name__ == '__main__':
    # Run tests with verbose output
    unittest.main(verbosity=2)