#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Copyright (C) 2013 University of Dundee & Open Microscopy Environment.
# All rights reserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#
#

import unittest
from ricerca import content



class TestRicerca(unittest.TestCase):
    def setUp(self):
        self.scale = 0.16125

    def createContentDB(self):
        cdb = {}
        cdb['info'] = {
            'name': 'omepslid2.compbio.cs.cmu.edu-slf0.pkl',
            'num_channels_required': 1,
            'slf_code': 'def calc_slf0(protimagef, scale):\n    return None\n',
            'slf_name': 'slf0'
            }
        cdb[self.scale] = [
            [
                1, 'host', 'user', 'metadata_url', 'img_url', 'render_url',
                1, 0, 0, 0, 0,
                5.1, 5.2,
                ],
            [
                2, 'host', 'user', 'metadata_url', 'img_url', 'render_url',
                2, 0, 0, 0, 0,
                1.0, 1.0,
                ],
            [
                3, 'host', 'user', 'metadata_url', 'img_url', 'render_url',
                3, 0, 0, 0, 0,
                1.2, 0.9,
                ],
            [
                4, 'host', 'user', 'metadata_url', 'img_url', 'render_url',
                4, 0, 0, 0, 0,
                1.1, 1.2,
                ],
            [
                5, 'host', 'user', 'metadata_url', 'img_url', 'render_url',
                5, 0, 0, 0, 0,
                8.1, 8.2,
                ],
            ]
        return cdb

    def createMiniTrainTest(self):
        train = [
            ['a', 0, [1., 2.]],
            ['b', 1, [5., 1.]],
            ]
        test = [
            ['c', 2, [-1., -2.]],
            ['d', 3, [-5., -1.]],
            ]
        return train, test

    def createCandidateGoodSet(self):
        candidates = [
            ['1.0.0.0.0', 0, [5.1, 5.2]],
            ['2.0.0.0.0', 0, [1.0, 1.0]],
            ]
        goodset = [
            ['3.0.0.0.0', 1, [1.2, 0.9]],
            ['4.0.0.0.0', 1, [1.1, 1.2]],
            ]
        return candidates, goodset

    def processIds(self, cdbr):
        return ['.'.join(str(c) for c in cdbr[6:11]), cdbr[2], cdbr[1]]

    def processSearchSet(self, cdb, im_ref_dict, dscale):
        self.assertIn(dscale, cdb)
        iid_cdb_dict = dict((k[6], k) for k in cdb[dscale])
        goodset_pos = []

        for id in im_ref_dict:
            iid = long(id.split('.')[0])
            self.assertIn(iid, iid_cdb_dict)
            feats = iid_cdb_dict[iid][11:]
            goodset_pos.append([id, 1, feats])
 
        #print 'im_ref_dict', im_ref_dict
        #print 'goodset_pos', goodset_pos

        return goodset_pos


    def test_getDBscales(self):
        db = { 'info': None, 1.5: [], 3.5: [] }
        s = content.getDBscales(db, 2.5)
        self.assertEqual(len(s), 1)
        self.assertEqual(s[0], 3.5)


    def test_rankingWrapper(self):
        cdb = self.createContentDB()
        scale = self.scale * 1.1
        im_ref_dict = {
            '2.0.0.0.0': [(scale, ''), 1],
            '5.0.0.0.0': [(scale, ''), -1]
            }

        #print 'processIds'
        #print '\n'.join(str(self.processIds(r)) for r in cdb[self.scale])

        #print 'processSearchSet'
        #print self.processSearchSet(cdb, im_ref_dict, self.scale)

        #print 'rankingWrapper'
        final_result, dscale = content.rankingWrapper(
            cdb, im_ref_dict, self.processIds, self.processSearchSet)
        #print 'final_result', final_result
        #print 'dscale', dscale

        self.assertEqual(dscale, self.scale)
        self.assertEqual([f[0] for f in final_result],
                         ['2.0.0.0.0',
                          '3.0.0.0.0',
                          '4.0.0.0.0',
                          '1.0.0.0.0',
                          '5.0.0.0.0'])




    def test_distance(self):
        alpha = .5
        candidates, goodset = self.createCandidateGoodSet()

        d = content.distance(alpha, candidates[0], goodset)
        self.assertAlmostEqual(d, 1.5472243)

        d = content.distance(alpha, candidates[1], goodset)
        self.assertAlmostEqual(d, 0.6876560)

    def test_norm(self):
        n = content.norm([2, 3], [5, 4], 2)
        # should equal sqrt(10)
        self.assertAlmostEqual(n, 3.1622777)

    @unittest.expectedFailure
    def test_featnorm(self):
        trainset, testset = self.createMiniTrainTest()

        # Not sure if this function is correct, since it can return values
        # <0 when it sounds like it should rescale to [0,1]
        # ... in which case if ther are only two values you'd expect them to
        # be 0 and 1 (or identical)
        train, test = content.featnorm(trainset, testset)
        self.assertEqual(train[0][:2], ['a', 0])
        self.assertEqual(train[1][:2], ['b', 1])
        self.assertEqual(test[0][:2], ['c', 2])
        self.assertEqual(test[1][:2], ['d', 3])

        def between01(xs):
            return all(x >= 0 and x <= 1 for x in xs)

        self.assertTrue(between01(train[0][2]))
        self.assertTrue(between01(train[1][2]))
        self.assertTrue(between01(test[0][2]))
        self.assertTrue(between01(test[1][2]))

        # TODO: figure out what this function should actually return


    def test_featnorm_z(self):
        trainset, testset = self.createMiniTrainTest()
        train, test = content.featnorm_z(trainset, testset)

        self.assertEqual(train[0][:2], ['a', 0])
        self.assertEqual(train[1][:2], ['b', 1])
        self.assertEqual(test[0][:2], ['c', 2])
        self.assertEqual(test[1][:2], ['d', 3])

        # TODO: figure out what this function should actually return


    def test_ranking(self):
        alpha = .5
        candidates, goodset = self.createCandidateGoodSet()

        normalization = ''
        r = content.ranking(alpha, candidates, goodset, normalization)
        self.assertEqual(r[0], ['2.0.0.0.0', '1.0.0.0.0'])
        self.assertAlmostEqual(r[1][0], 0.4564908)
        self.assertAlmostEqual(r[1][1], 1.0269875)

        normalization = 'zscore'
        r = content.ranking(alpha, candidates, goodset, normalization)
        self.assertEqual(r[0], ['2.0.0.0.0', '1.0.0.0.0'])
        self.assertAlmostEqual(r[1][0], 0.5729756)
        self.assertAlmostEqual(r[1][1], 1.2890159)

        # TODO: Check the numbers are actually correct




if __name__ == '__main__':
    unittest.main()
