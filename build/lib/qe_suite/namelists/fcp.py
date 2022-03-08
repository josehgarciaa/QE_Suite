from .namelist import Namelist


class FCP(Namelist):

    def __init__(self):
        self.set_parameters(
            {"fcp_mu": None, "fcp_dynamics": None, "fcp_conv_thr": None,
             "fcp_ndiis": None, "fcp_mass": None, "fcp_velocity": None,
             "fcp_temperature": None, "fcp_tempw": None, "fcp_tolp": None,
             "fcp_delta_t": None, "fcp_nraise": None, "freeze_all_atoms": None}
        )
        self.set_name("&FCP");
