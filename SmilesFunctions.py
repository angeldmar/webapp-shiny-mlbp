from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors3D
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import Draw
import pandas as pd


def descriptor_string_partitioner(dataframe_raw, selected_descriptor):

    dataframe_subset = dataframe_raw.astype('string')
    dataframe_subset[selected_descriptor] = dataframe_subset[selected_descriptor].str.replace("[", "", regex=True).str.replace("]", "", regex=True)
    dataframe_subset = dataframe_subset[selected_descriptor].str.split(pat=",", expand=True).add_prefix(selected_descriptor + '_value_').fillna(0)
    final_dataframe = dataframe_raw.join(dataframe_subset)
    final_dataframe = final_dataframe.drop(selected_descriptor, axis=1)

    return final_dataframe


def dictionary_to_dataframe(dictionary_smiles):

    temporary_array_maker = {"molecule": dictionary_smiles}
    smiles_dataframe = pd.DataFrame.from_dict(temporary_array_maker)
    smiles_dataframe = smiles_dataframe.transpose()

    return smiles_dataframe

def descriptor3d_data_transformation(descriptors3d_raw):

    descriptors3d_dataframe = dictionary_to_dataframe(descriptors3d_raw)
    final_descriptors3d_dataframe = descriptor_string_partitioner(descriptors3d_dataframe, 'CalcAUTOCORR3D')
    final_descriptors3d_dataframe = descriptor_string_partitioner(final_descriptors3d_dataframe, 'CalcRDF')
    final_descriptors3d_dataframe = descriptor_string_partitioner(final_descriptors3d_dataframe, 'CalcMORSE')
    final_descriptors3d_dataframe = descriptor_string_partitioner(final_descriptors3d_dataframe, 'CalcWHIM')
    final_descriptors3d_dataframe = descriptor_string_partitioner(final_descriptors3d_dataframe, 'CalcGETAWAY')

    return final_descriptors3d_dataframe

def descriptor2d_data_transformation(descriptors2d_raw):

    descriptors2d_dataframe = dictionary_to_dataframe(descriptors2d_raw)
    final_descriptors2d_dataframe = descriptor_string_partitioner(descriptors2d_dataframe, 'MQNs_')
    final_descriptors2d_dataframe = descriptor_string_partitioner(final_descriptors2d_dataframe, 'CalcAUTOCORR2D')

    return final_descriptors2d_dataframe


def calculating_descriptors3d(string_smiles):

    descriptors3d_dictionary = {}
    m_in_function = Chem.MolFromSmiles(string_smiles)
    m2 = Chem.AddHs(m_in_function)
    AllChem.EmbedMolecule(m2)
    AllChem.MMFFOptimizeMolecule(m2)

    descriptors3d_dictionary['CalcPBF'] = rdMolDescriptors.CalcPBF(m2)
    descriptors3d_dictionary['CalcAUTOCORR3D'] = rdMolDescriptors.CalcAUTOCORR3D(m2)
    descriptors3d_dictionary['CalcRDF'] = rdMolDescriptors.CalcRDF(m2)
    descriptors3d_dictionary['CalcMORSE'] = rdMolDescriptors.CalcMORSE(m2)
    descriptors3d_dictionary['CalcWHIM'] = rdMolDescriptors.CalcWHIM(m2)
    descriptors3d_dictionary['CalcGETAWAY'] = rdMolDescriptors.CalcGETAWAY(m2)
    descriptors3d_dictionary['Asphericity'] = Descriptors3D.Asphericity(m2)
    descriptors3d_dictionary['Eccentricity'] = Descriptors3D.Eccentricity(m2)
    descriptors3d_dictionary['InertialShapeFactor'] = Descriptors3D.InertialShapeFactor(m2)
    descriptors3d_dictionary['NPR1'] = Descriptors3D.NPR1(m2)
    descriptors3d_dictionary['NPR2'] = Descriptors3D.NPR2(m2)
    descriptors3d_dictionary['PMI1'] = Descriptors3D.PMI1(m2)
    descriptors3d_dictionary['PMI2'] = Descriptors3D.PMI2(m2)
    descriptors3d_dictionary['PMI3'] = Descriptors3D.PMI3(m2)
    descriptors3d_dictionary['RadiusOfGyration'] = Descriptors3D.RadiusOfGyration(m2)
    descriptors3d_dictionary['SpherocityIndex'] = Descriptors3D.SpherocityIndex(m2)

    return descriptors3d_dictionary


def calculating_descriptors2d(string_smiles):

    descriptors2d_dictionary = {}
    m = Chem.MolFromSmiles(string_smiles)

    descriptors2d_dictionary['BalabanJ'] = Descriptors.BalabanJ(m)
    descriptors2d_dictionary['BertzCT'] = Descriptors.BertzCT(m)
    descriptors2d_dictionary['Ipc'] = Descriptors.Ipc(m)
    descriptors2d_dictionary['HallKierAlpha'] = Descriptors.HallKierAlpha(m)
    descriptors2d_dictionary['Kappa1'] = Descriptors.Kappa1(m)
    descriptors2d_dictionary['Kappa2'] = Descriptors.Kappa2(m)
    descriptors2d_dictionary['Kappa3'] = Descriptors.Kappa3(m)
    descriptors2d_dictionary['Chi0'] = Descriptors.Chi0(m)
    descriptors2d_dictionary['Chi1'] = Descriptors.Chi1(m)
    descriptors2d_dictionary['Chi0n'] = Descriptors.Chi0n(m)
    descriptors2d_dictionary['Chi1n'] = Descriptors.Chi1n(m)
    descriptors2d_dictionary['Chi2n'] = Descriptors.Chi2n(m)
    descriptors2d_dictionary['Chi3n'] = Descriptors.Chi3n(m)
    descriptors2d_dictionary['Chi0v'] = Descriptors.Chi0v(m)
    descriptors2d_dictionary['Chi1v'] = Descriptors.Chi1v(m)
    descriptors2d_dictionary['Chi2v'] = Descriptors.Chi2v(m)
    descriptors2d_dictionary['Chi3v'] = Descriptors.Chi3v(m)
    descriptors2d_dictionary['Chi4v'] = Descriptors.Chi4v(m)
    descriptors2d_dictionary['MolLogP'] = Descriptors.MolLogP(m)
    descriptors2d_dictionary['MolMR'] = Descriptors.MolMR(m)
    descriptors2d_dictionary['MolWt'] = Descriptors.MolWt(m)
    descriptors2d_dictionary['ExactMolWt'] = Descriptors.ExactMolWt(m)
    descriptors2d_dictionary['HeavyAtomCount'] = Descriptors.HeavyAtomCount(m)
    descriptors2d_dictionary['HeavyAtomMolWt'] = Descriptors.HeavyAtomMolWt(m)
    descriptors2d_dictionary['NHOHCount'] = Descriptors.NHOHCount(m)
    descriptors2d_dictionary['NOCount'] = Descriptors.NOCount(m)
    descriptors2d_dictionary['NumHAcceptors'] = Descriptors.NumHAcceptors(m)
    descriptors2d_dictionary['NumHDonors'] = Descriptors.NumHDonors(m)
    descriptors2d_dictionary['NumHeteroatoms'] = Descriptors.NumHeteroatoms(m)
    descriptors2d_dictionary['NumRotableBonds'] = Descriptors.NumRotatableBonds(m)
    descriptors2d_dictionary['NumValenceElectrons'] = Descriptors.NumValenceElectrons(m)
    descriptors2d_dictionary['NumAromaticRings'] = Descriptors.NumAromaticRings(m)
    descriptors2d_dictionary['NumSaturatedRings'] = Descriptors.NumSaturatedRings(m)
    descriptors2d_dictionary['NumAliphaticRings'] = Descriptors.NumAliphaticRings(m)
    descriptors2d_dictionary['NumAromaticHeterocycles'] = Descriptors.NumAromaticHeterocycles(m)
    descriptors2d_dictionary['NumSaturatedHeterocycles'] = Descriptors.NumSaturatedHeterocycles(m)
    descriptors2d_dictionary['NumSaturatedCarbocycles'] = Descriptors.NumSaturatedCarbocycles(m)
    descriptors2d_dictionary['NumAliphaticHeterocyles'] = Descriptors.NumAliphaticHeterocycles(m)
    descriptors2d_dictionary['NumAliphaticCarbocycles'] = Descriptors.NumAliphaticCarbocycles(m)
    descriptors2d_dictionary['RingCount'] = Descriptors.RingCount(m)
    descriptors2d_dictionary['FractionCSP3'] = Descriptors.FractionCSP3(m)
    descriptors2d_dictionary['TPSA'] = Descriptors.TPSA(m)
    descriptors2d_dictionary['LabuteASA'] = Descriptors.LabuteASA(m)
    descriptors2d_dictionary['PEOE_VSA1'] = Descriptors.PEOE_VSA1(m)
    descriptors2d_dictionary['PEOE_VSA2'] = Descriptors.PEOE_VSA2(m)
    descriptors2d_dictionary['PEOE_VSA3'] = Descriptors.PEOE_VSA3(m)
    descriptors2d_dictionary['PEOE_VSA4'] = Descriptors.PEOE_VSA4(m)
    descriptors2d_dictionary['PEOE_VSA5'] = Descriptors.PEOE_VSA5(m)
    descriptors2d_dictionary['PEOE_VSA6'] = Descriptors.PEOE_VSA6(m)
    descriptors2d_dictionary['PEOE_VSA7'] = Descriptors.PEOE_VSA7(m)
    descriptors2d_dictionary['PEOE_VSA8'] = Descriptors.PEOE_VSA8(m)
    descriptors2d_dictionary['PEOE_VSA9'] = Descriptors.PEOE_VSA9(m)
    descriptors2d_dictionary['PEOE_VSA10'] = Descriptors.PEOE_VSA10(m)
    descriptors2d_dictionary['PEOE_VSA11'] = Descriptors.PEOE_VSA11(m)
    descriptors2d_dictionary['PEOE_VSA12'] = Descriptors.PEOE_VSA12(m)
    descriptors2d_dictionary['PEOE_VSA13'] = Descriptors.PEOE_VSA13(m)
    descriptors2d_dictionary['PEOE_VSA14'] = Descriptors.PEOE_VSA14(m)
    descriptors2d_dictionary['SMR_VSA1'] = Descriptors.SMR_VSA1(m)
    descriptors2d_dictionary['SMR_VSA2'] = Descriptors.SMR_VSA2(m)
    descriptors2d_dictionary['SMR_VSA3'] = Descriptors.SMR_VSA3(m)
    descriptors2d_dictionary['SMR_VSA4'] = Descriptors.SMR_VSA4(m)
    descriptors2d_dictionary['SMR_VSA5'] = Descriptors.SMR_VSA5(m)
    descriptors2d_dictionary['SMR_VSA6'] = Descriptors.SMR_VSA6(m)
    descriptors2d_dictionary['SMR_VSA7'] = Descriptors.SMR_VSA7(m)
    descriptors2d_dictionary['SMR_VSA8'] = Descriptors.SMR_VSA8(m)
    descriptors2d_dictionary['SMR_VSA9'] = Descriptors.SMR_VSA9(m)
    descriptors2d_dictionary['SMR_VSA10'] = Descriptors.SMR_VSA10(m)
    descriptors2d_dictionary['SlogP_VSA1'] = Descriptors.SlogP_VSA1(m)
    descriptors2d_dictionary['SlogP_VSA2'] = Descriptors.SlogP_VSA2(m)
    descriptors2d_dictionary['SlogP_VSA3'] = Descriptors.SlogP_VSA3(m)
    descriptors2d_dictionary['SlogP_VSA4'] = Descriptors.SlogP_VSA4(m)
    descriptors2d_dictionary['SlogP_VSA5'] = Descriptors.SlogP_VSA5(m)
    descriptors2d_dictionary['SlogP_VSA6'] = Descriptors.SlogP_VSA6(m)
    descriptors2d_dictionary['SlogP_VSA7'] = Descriptors.SlogP_VSA7(m)
    descriptors2d_dictionary['SlogP_VSA8'] = Descriptors.SlogP_VSA8(m)
    descriptors2d_dictionary['SlogP_VSA9'] = Descriptors.SlogP_VSA9(m)
    descriptors2d_dictionary['SlogP_VSA10'] = Descriptors.SlogP_VSA10(m)
    descriptors2d_dictionary['SlogP_VSA11'] = Descriptors.SlogP_VSA11(m)
    descriptors2d_dictionary['SlogP_VSA12'] = Descriptors.SlogP_VSA12(m)
    descriptors2d_dictionary['EState_VSA1'] = Descriptors.EState_VSA1(m)
    descriptors2d_dictionary['EState_VSA2'] = Descriptors.EState_VSA2(m)
    descriptors2d_dictionary['EState_VSA3'] = Descriptors.EState_VSA3(m)
    descriptors2d_dictionary['EState_VSA4'] = Descriptors.EState_VSA4(m)
    descriptors2d_dictionary['EState_VSA5'] = Descriptors.EState_VSA5(m)
    descriptors2d_dictionary['EState_VSA6'] = Descriptors.EState_VSA6(m)
    descriptors2d_dictionary['EState_VSA7'] = Descriptors.EState_VSA7(m)
    descriptors2d_dictionary['EState_VSA8'] = Descriptors.EState_VSA8(m)
    descriptors2d_dictionary['EState_VSA9'] = Descriptors.EState_VSA9(m)
    descriptors2d_dictionary['EState_VSA10'] = Descriptors.EState_VSA10(m)
    descriptors2d_dictionary['EState_VSA11'] = Descriptors.EState_VSA11(m)
    descriptors2d_dictionary['VSA_EState1'] = Descriptors.VSA_EState1(m)
    descriptors2d_dictionary['VSA_EState2'] = Descriptors.VSA_EState2(m)
    descriptors2d_dictionary['VSA_EState3'] = Descriptors.VSA_EState3(m)
    descriptors2d_dictionary['VSA_EState4'] = Descriptors.VSA_EState4(m)
    descriptors2d_dictionary['VSA_EState5'] = Descriptors.VSA_EState5(m)
    descriptors2d_dictionary['VSA_EState6'] = Descriptors.VSA_EState6(m)
    descriptors2d_dictionary['VSA_EState7'] = Descriptors.VSA_EState7(m)
    descriptors2d_dictionary['VSA_EState8'] = Descriptors.VSA_EState8(m)
    descriptors2d_dictionary['VSA_EState9'] = Descriptors.VSA_EState9(m)
    descriptors2d_dictionary['VSA_EState10'] = Descriptors.VSA_EState10(m)
    descriptors2d_dictionary['CalcPhi'] = rdMolDescriptors.CalcPhi(m)
    descriptors2d_dictionary['NumAmideBonds'] = rdMolDescriptors.CalcNumAmideBonds(m)
    descriptors2d_dictionary['Chi4n'] = rdMolDescriptors.CalcChi4n(m)
    descriptors2d_dictionary['CalcNumAromaticCarbocycles'] = rdMolDescriptors.CalcNumAromaticCarbocycles(m)
    descriptors2d_dictionary['CalcNumSpiroAtoms'] = rdMolDescriptors.CalcNumSpiroAtoms(m)
    descriptors2d_dictionary['CalcNumBridgeheadAtoms'] = rdMolDescriptors.CalcNumBridgeheadAtoms(m)
    descriptors2d_dictionary['MQNs_'] = rdMolDescriptors.MQNs_(m)
    descriptors2d_dictionary['CalcAUTOCORR2D'] = rdMolDescriptors.CalcAUTOCORR2D(m)

    return descriptors2d_dictionary

def calculating_lipinski_dataframe(string_smiles):

    lipinski_dictionary = {}
    m = Chem.MolFromSmiles(string_smiles)

    lipinski_dictionary['FractionCSP3'] = Lipinski.FractionCSP3(m)
    lipinski_dictionary['HeavyAtomCount'] = Lipinski.HeavyAtomCount(m)
    lipinski_dictionary['NHOHCount'] = Lipinski.NHOHCount(m)
    lipinski_dictionary['NOCount'] = Lipinski.NOCount(m)
    lipinski_dictionary['NumAliphaticCarbocycles'] = Lipinski.NumAliphaticCarbocycles(m)
    lipinski_dictionary['NumAliphaticHeterocycles'] = Lipinski.NumAliphaticHeterocycles(m)
    lipinski_dictionary['NumAromaticCarbocycles'] = Lipinski.NumAromaticCarbocycles(m)
    lipinski_dictionary['NumHAcceptors'] = Lipinski.NumHAcceptors(m)
    lipinski_dictionary['NumHDonors'] = Lipinski.NumHDonors(m)
    lipinski_dictionary['NumHeteroatoms'] = Lipinski.NumHeteroatoms(m)
    lipinski_dictionary['NumRotatableBonds'] = Lipinski.NumRotatableBonds(m)
    lipinski_dictionary['NumSaturatedCarbocycles'] = Lipinski.NumSaturatedCarbocycles(m)
    lipinski_dictionary['NumSaturatedHeterocycles'] = Lipinski.NumSaturatedHeterocycles(m)
    lipinski_dictionary['NumSaturatedRings'] = Lipinski.NumSaturatedRings(m)
    lipinski_dictionary['RingCount'] = Lipinski.RingCount(m)

    lipinski_dataframe = dictionary_to_dataframe(lipinski_dictionary)

    return lipinski_dataframe


def drawing_smiles(string_smiles):
    m = Chem.MolFromSmiles(string_smiles)
    AllChem.Compute2DCoords(m)
    Draw.MolToFile(m, "2D_smiles.o.png", size=(1280,720))


#  Entry point Solo debe usarse en test, de lo contrario generara errores en reticulate


if __name__ == "__main__":
    # smiles_de_prueba = input("Introduce tu smiles:")
    # descriptor3d_de_prueba = calculating_descriptors3d(smiles_de_prueba)
    # descriptor3d_de_prueba_transformado = descriptor3d_data_transformation(descriptor3d_de_prueba)
    # descriptor2d_de_prueba = calculating_descriptors2d(smiles_de_prueba)
    # descriptor2d_de_prueba_transformado = descriptor2d_data_transformation(descriptor2d_de_prueba)
    # descriptores_lipinski = calculating_lipinski_dataframe(smiles_de_prueba)
    # print(descriptor3d_de_prueba)
    # print(descriptor3d_de_prueba_transformado)
    # print(descriptor2d_de_prueba)
    # print(descriptor2d_de_prueba_transformado)
    # print(descriptores_lipinski)
    drawing_smiles("CCC[C@@H](O)CC\C=C\C=C\C#CC#C\C=C\CO")
