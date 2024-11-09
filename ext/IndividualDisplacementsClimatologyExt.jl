module IndividualDisplacementsClimatologyExt

using Climatology, IndividualDisplacements
import IndividualDisplacements: data_path


function data_path(name::Symbol)
    if name==:OCCA
        Climatology.get_occa_velocity_if_needed()
        Climatology.ScratchSpaces.OCCA
    elseif name==:ECCO
        Climatology.get_ecco_velocity_if_needed()
        Climatology.get_ecco_variable_if_needed("THETA")
        Climatology.get_ecco_variable_if_needed("SALT")
        Climatology.ScratchSpaces.ECCO        
    else
        tempdir()
    end
end

end
