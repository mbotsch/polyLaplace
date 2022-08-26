#include "InterfaceAttrib.hh"
#include "../Core/TopologyKernel.hh"

namespace OpenVolumeMesh {

InterfaceAttrib::InterfaceAttrib(TopologyKernel &_kernel)
    : kernel_(_kernel),
      v_interface_(&_kernel, "interface", false),
      e_interface_(&_kernel, "interface", false),
      f_interface_(&_kernel, "interface", false)
{
    _kernel.set_persistent(v_interface_, true);
    _kernel.set_persistent(e_interface_, true);
    _kernel.set_persistent(f_interface_, true);
}

}
