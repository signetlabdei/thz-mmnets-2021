# -*- Mode: python; py-indent-offset: 4; indent-tabs-mode: nil; coding: utf-8; -*-

# def options(opt):
#     pass

# def configure(conf):
#     conf.check_nonfatal(header_name='stdint.h', define_name='HAVE_STDINT_H')

def build(bld):
    module = bld.create_ns3_module('thz', ['core', 'propagation', 'internet', 'spectrum', 'applications', 'mobility', 'antenna', 'network','qd-channel'])
    module.source = [
               
        'model/thz-spectrum-propagation-loss.cc',
        'model/thz-channel.cc',
        'model/thz-dir-antenna.cc',
        'model/thz-spectrum-waveform.cc',
        'model/thz-spectrum-signal-parameters.cc',
        'model/thz-phy-nano.cc',
        'model/thz-phy-macro.cc',
        'model/thz-mac-macro.cc',
        'model/thz-mac-macro-ap.cc',
        'model/thz-mac-macro-client.cc',
        'model/thz-mac-nano.cc',
        'model/thz-mac-header.cc',
        'model/thz-net-device.cc',
        'model/thz-energy-model.cc',
        'model/traffic-generator.cc',
        'model/thz-udp-server.cc',
        'model/thz-udp-client.cc',
        'model/thz-udp-trace-client.cc',
                     
        'helper/thz-helper.cc',
        'helper/thz-mac-nano-helper.cc',
        'helper/thz-mac-macro-helper.cc',
        'helper/thz-mac-macro-ap-helper.cc',
        'helper/thz-mac-macro-client-helper.cc',
        'helper/thz-phy-nano-helper.cc',
        'helper/thz-phy-macro-helper.cc',
        'helper/thz-directional-antenna-helper.cc',
        'helper/traffic-generator-helper.cc',
        'helper/thz-energy-model-helper.cc',
        'helper/thz-udp-client-server-helper.cc',
        ]

    module_test = bld.create_ns3_module_test_library('thz')
    module_test.source = [
        'test/test-thz-directional-antenna.cc',
        'test/test-thz-path-loss.cc',
        'test/test-thz-psd-macro.cc',
        'test/test-thz-psd-nano.cc',
        ]

    headers = bld(features='ns3header')
    headers.module = 'thz'
    headers.source = [
        'model/thz-spectrum-propagation-loss.h',
        'model/thz-channel.h',
        'model/thz-dir-antenna.h',
        'model/thz-spectrum-waveform.h',
        'model/thz-spectrum-signal-parameters.h',
        'model/thz-phy.h',
        'model/thz-phy-nano.h',
        'model/thz-phy-macro.h',
        'model/thz-mac.h',
        'model/thz-mac-macro.h',
        'model/thz-mac-macro-ap.h',
        'model/thz-mac-macro-client.h',
        'model/thz-mac-nano.h',
        'model/thz-mac-header.h',
        'model/thz-net-device.h',
        'model/thz-energy-model.h',
        'model/traffic-generator.h',
        'model/thz-udp-server.h',
        'model/thz-udp-client.h',
        'model/thz-udp-trace-client.h',

        'helper/thz-helper.h',
        'helper/thz-mac-nano-helper.h',
        'helper/thz-mac-macro-helper.h',
        'helper/thz-mac-macro-ap-helper.h',
        'helper/thz-mac-macro-client-helper.h',
        'helper/thz-phy-nano-helper.h',
        'helper/thz-phy-macro-helper.h',
        'helper/thz-directional-antenna-helper.h',
        'helper/traffic-generator-helper.h',
        'helper/thz-energy-model-helper.h',
        'helper/thz-udp-client-server-helper.h',
        ]

    if bld.env.ENABLE_EXAMPLES:
        bld.recurse('examples')

    # bld.ns3_python_bindings()

