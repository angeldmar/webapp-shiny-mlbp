from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors3D
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Lipinski
import pandas as pd


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


def descriptor_string_partitioner(dataframe_raw, selected_descriptor):
    descriptors3d_dataframe_subset = dataframe_raw.astype('string')
    descriptors3d_dataframe_subset[selected_descriptor] = descriptors3d_dataframe_subset[selected_descriptor].str.replace("[", "").str.replace("]", "")
    descriptors3d_dataframe_subset = descriptors3d_dataframe_subset[selected_descriptor].str.split(pat=",", expand=True).add_prefix(selected_descriptor + '_value_').fillna(0)
    final_dataframe = dataframe_raw.join(descriptors3d_dataframe_subset)
    final_dataframe = final_dataframe.drop(selected_descriptor, axis=1)
    return final_dataframe


def descriptor3d_data_transformation(descriptors3d_raw):
    temporary_array_maker = {"molecule": descriptors3d_raw}
    descriptors3d_dataframe = pd.DataFrame.from_dict(temporary_array_maker)
    descriptors3d_dataframe = descriptors3d_dataframe.transpose()

    final_descriptors3d_dataframe = descriptor_string_partitioner(descriptors3d_dataframe, 'CalcAUTOCORR3D')
    final_descriptors3d_dataframe = descriptor_string_partitioner(descriptors3d_dataframe, 'CalcRDF')
    final_descriptors3d_dataframe = descriptor_string_partitioner(descriptors3d_dataframe, 'CalcMORSE')
    final_descriptors3d_dataframe = descriptor_string_partitioner(descriptors3d_dataframe, 'CalcWHIM')
    final_descriptors3d_dataframe = descriptor_string_partitioner(descriptors3d_dataframe, 'CalcGETAWAY')

    return final_descriptors3d_dataframe

# Descriptores3D_Dataframe.to_csv("./DescriptorsCSV/Descriptores3D_Dataframe.csv")




    # Exportando base de datos de smiles a una lista
def run():
    smiles = pd.read_csv('./SMILES.csv', delimiter=',')
    smiles = smiles.transpose()
    smiles = smiles.values.tolist()

    # Recorriendo la lista para calcular descriptores

    Descriptores2D = {}
    for x in smiles:
        for molecula in x:
            Descriptores2D_diccionario = {}
            m = Chem.MolFromSmiles(molecula)

            Descriptores2D_diccionario['BalabanJ'] = Descriptors.BalabanJ(m)
            Descriptores2D_diccionario['BertzCT'] = Descriptors.BertzCT(m)
            Descriptores2D_diccionario['Ipc'] = Descriptors.Ipc(m)
            Descriptores2D_diccionario['HallKierAlpha'] = Descriptors.HallKierAlpha(m)
            Descriptores2D_diccionario['Kappa1'] = Descriptors.Kappa1(m)
            Descriptores2D_diccionario['Kappa2'] = Descriptors.Kappa2(m)
            Descriptores2D_diccionario['Kappa3'] = Descriptors.Kappa3(m)
            Descriptores2D_diccionario['Chi0'] = Descriptors.Chi0(m)
            Descriptores2D_diccionario['Chi1'] = Descriptors.Chi1(m)
            Descriptores2D_diccionario['Chi0n'] = Descriptors.Chi0n(m)
            Descriptores2D_diccionario['Chi1n'] = Descriptors.Chi1n(m)
            Descriptores2D_diccionario['Chi2n'] = Descriptors.Chi2n(m)
            Descriptores2D_diccionario['Chi3n'] = Descriptors.Chi3n(m)
            Descriptores2D_diccionario['Chi0v'] = Descriptors.Chi0v(m)
            Descriptores2D_diccionario['Chi1v'] = Descriptors.Chi1v(m)
            Descriptores2D_diccionario['Chi2v'] = Descriptors.Chi2v(m)
            Descriptores2D_diccionario['Chi3v'] = Descriptors.Chi3v(m)
            Descriptores2D_diccionario['Chi4v'] = Descriptors.Chi4v(m)
            Descriptores2D_diccionario['MolLogP'] = Descriptors.MolLogP(m)
            Descriptores2D_diccionario['MolMR'] = Descriptors.MolMR(m)
            Descriptores2D_diccionario['MolWt'] = Descriptors.MolWt(m)
            Descriptores2D_diccionario['ExactMolWt'] = Descriptors.ExactMolWt(m)
            Descriptores2D_diccionario['HeavyAtomCount'] = Descriptors.HeavyAtomCount(m)
            Descriptores2D_diccionario['HeavyAtomMolWt'] = Descriptors.HeavyAtomMolWt(m)
            Descriptores2D_diccionario['NHOHCount'] = Descriptors.NHOHCount(m)
            Descriptores2D_diccionario['NOCount'] = Descriptors.NOCount(m)
            Descriptores2D_diccionario['NumHAcceptors'] = Descriptors.NumHAcceptors(m)
            Descriptores2D_diccionario['NumHDonors'] = Descriptors.NumHDonors(m)
            Descriptores2D_diccionario['NumHeteroatoms'] = Descriptors.NumHeteroatoms(m)
            Descriptores2D_diccionario['NumRotableBonds'] = Descriptors.NumRotatableBonds(m)
            Descriptores2D_diccionario['NumValenceElectrons'] = Descriptors.NumValenceElectrons(m)
            Descriptores2D_diccionario['NumAromaticRings'] = Descriptors.NumAromaticRings(m)
            Descriptores2D_diccionario['NumSaturatedRings'] = Descriptors.NumSaturatedRings(m)
            Descriptores2D_diccionario['NumAliphaticRings'] = Descriptors.NumAliphaticRings(m)
            Descriptores2D_diccionario['NumAromaticHeterocycles'] = Descriptors.NumAromaticHeterocycles(m)
            Descriptores2D_diccionario['NumSaturatedHeterocycles'] = Descriptors.NumSaturatedHeterocycles(m)
            Descriptores2D_diccionario['NumSaturatedCarbocycles'] = Descriptors.NumSaturatedCarbocycles(m)
            Descriptores2D_diccionario['NumAliphaticHeterocyles'] = Descriptors.NumAliphaticHeterocycles(m)
            Descriptores2D_diccionario['NumAliphaticCarbocycles'] = Descriptors.NumAliphaticCarbocycles(m)
            Descriptores2D_diccionario['RingCount'] = Descriptors.RingCount(m)
            Descriptores2D_diccionario['FractionCSP3'] = Descriptors.FractionCSP3(m)
            Descriptores2D_diccionario['TPSA'] = Descriptors.TPSA(m)
            Descriptores2D_diccionario['LabuteASA'] = Descriptors.LabuteASA(m)
            Descriptores2D_diccionario['PEOE_VSA1'] = Descriptors.PEOE_VSA1(m)
            Descriptores2D_diccionario['PEOE_VSA2'] = Descriptors.PEOE_VSA2(m)
            Descriptores2D_diccionario['PEOE_VSA3'] = Descriptors.PEOE_VSA3(m)
            Descriptores2D_diccionario['PEOE_VSA4'] = Descriptors.PEOE_VSA4(m)
            Descriptores2D_diccionario['PEOE_VSA5'] = Descriptors.PEOE_VSA5(m)
            Descriptores2D_diccionario['PEOE_VSA6'] = Descriptors.PEOE_VSA6(m)
            Descriptores2D_diccionario['PEOE_VSA7'] = Descriptors.PEOE_VSA7(m)
            Descriptores2D_diccionario['PEOE_VSA8'] = Descriptors.PEOE_VSA8(m)
            Descriptores2D_diccionario['PEOE_VSA9'] = Descriptors.PEOE_VSA9(m)
            Descriptores2D_diccionario['PEOE_VSA10'] = Descriptors.PEOE_VSA10(m)
            Descriptores2D_diccionario['PEOE_VSA11'] = Descriptors.PEOE_VSA11(m)
            Descriptores2D_diccionario['PEOE_VSA12'] = Descriptors.PEOE_VSA12(m)
            Descriptores2D_diccionario['PEOE_VSA13'] = Descriptors.PEOE_VSA13(m)
            Descriptores2D_diccionario['PEOE_VSA14'] = Descriptors.PEOE_VSA14(m)
            Descriptores2D_diccionario['SMR_VSA1'] = Descriptors.SMR_VSA1(m)
            Descriptores2D_diccionario['SMR_VSA2'] = Descriptors.SMR_VSA2(m)
            Descriptores2D_diccionario['SMR_VSA3'] = Descriptors.SMR_VSA3(m)
            Descriptores2D_diccionario['SMR_VSA4'] = Descriptors.SMR_VSA4(m)
            Descriptores2D_diccionario['SMR_VSA5'] = Descriptors.SMR_VSA5(m)
            Descriptores2D_diccionario['SMR_VSA6'] = Descriptors.SMR_VSA6(m)
            Descriptores2D_diccionario['SMR_VSA7'] = Descriptors.SMR_VSA7(m)
            Descriptores2D_diccionario['SMR_VSA8'] = Descriptors.SMR_VSA8(m)
            Descriptores2D_diccionario['SMR_VSA9'] = Descriptors.SMR_VSA9(m)
            Descriptores2D_diccionario['SMR_VSA10'] = Descriptors.SMR_VSA10(m)
            Descriptores2D_diccionario['SlogP_VSA1'] = Descriptors.SlogP_VSA1(m)
            Descriptores2D_diccionario['SlogP_VSA2'] = Descriptors.SlogP_VSA2(m)
            Descriptores2D_diccionario['SlogP_VSA3'] = Descriptors.SlogP_VSA3(m)
            Descriptores2D_diccionario['SlogP_VSA4'] = Descriptors.SlogP_VSA4(m)
            Descriptores2D_diccionario['SlogP_VSA5'] = Descriptors.SlogP_VSA5(m)
            Descriptores2D_diccionario['SlogP_VSA6'] = Descriptors.SlogP_VSA6(m)
            Descriptores2D_diccionario['SlogP_VSA7'] = Descriptors.SlogP_VSA7(m)
            Descriptores2D_diccionario['SlogP_VSA8'] = Descriptors.SlogP_VSA8(m)
            Descriptores2D_diccionario['SlogP_VSA9'] = Descriptors.SlogP_VSA9(m)
            Descriptores2D_diccionario['SlogP_VSA10'] = Descriptors.SlogP_VSA10(m)
            Descriptores2D_diccionario['SlogP_VSA11'] = Descriptors.SlogP_VSA11(m)
            Descriptores2D_diccionario['SlogP_VSA12'] = Descriptors.SlogP_VSA12(m)
            Descriptores2D_diccionario['EState_VSA1'] = Descriptors.EState_VSA1(m)
            Descriptores2D_diccionario['EState_VSA2'] = Descriptors.EState_VSA2(m)
            Descriptores2D_diccionario['EState_VSA3'] = Descriptors.EState_VSA3(m)
            Descriptores2D_diccionario['EState_VSA4'] = Descriptors.EState_VSA4(m)
            Descriptores2D_diccionario['EState_VSA5'] = Descriptors.EState_VSA5(m)
            Descriptores2D_diccionario['EState_VSA6'] = Descriptors.EState_VSA6(m)
            Descriptores2D_diccionario['EState_VSA7'] = Descriptors.EState_VSA7(m)
            Descriptores2D_diccionario['EState_VSA8'] = Descriptors.EState_VSA8(m)
            Descriptores2D_diccionario['EState_VSA9'] = Descriptors.EState_VSA9(m)
            Descriptores2D_diccionario['EState_VSA10'] = Descriptors.EState_VSA10(m)
            Descriptores2D_diccionario['EState_VSA11'] = Descriptors.EState_VSA11(m)
            Descriptores2D_diccionario['VSA_EState1'] = Descriptors.VSA_EState1(m)
            Descriptores2D_diccionario['VSA_EState2'] = Descriptors.VSA_EState2(m)
            Descriptores2D_diccionario['VSA_EState3'] = Descriptors.VSA_EState3(m)
            Descriptores2D_diccionario['VSA_EState4'] = Descriptors.VSA_EState4(m)
            Descriptores2D_diccionario['VSA_EState5'] = Descriptors.VSA_EState5(m)
            Descriptores2D_diccionario['VSA_EState6'] = Descriptors.VSA_EState6(m)
            Descriptores2D_diccionario['VSA_EState7'] = Descriptors.VSA_EState7(m)
            Descriptores2D_diccionario['VSA_EState8'] = Descriptors.VSA_EState8(m)
            Descriptores2D_diccionario['VSA_EState9'] = Descriptors.VSA_EState9(m)
            Descriptores2D_diccionario['VSA_EState10'] = Descriptors.VSA_EState10(m)
            Descriptores2D_diccionario['CalcPhi'] = rdMolDescriptors.CalcPhi(m)
            Descriptores2D_diccionario['NumAmideBonds'] = rdMolDescriptors.CalcNumAmideBonds(m)
            Descriptores2D_diccionario['Chi4n'] = rdMolDescriptors.CalcChi4n(m)
            Descriptores2D_diccionario['CalcNumAromaticCarbocycles'] = rdMolDescriptors.CalcNumAromaticCarbocycles(m)
            Descriptores2D_diccionario['CalcNumSpiroAtoms'] = rdMolDescriptors.CalcNumSpiroAtoms(m)
            Descriptores2D_diccionario['CalcNumBridgeheadAtoms'] = rdMolDescriptors.CalcNumBridgeheadAtoms(m)
            Descriptores2D_diccionario['MQNs_'] = rdMolDescriptors.MQNs_(m)
            Descriptores2D_diccionario['CalcAUTOCORR2D'] = rdMolDescriptors.CalcAUTOCORR2D(m)

            Descriptores2D[molecula] = Descriptores2D_diccionario

    # Los descriptores almacenados en un diccionario se convierten a dataframes
    # Se realizan correcciones en el formato del dataframe
    # Se guardan como csv

    Descriptores2D_Dataframe = pd.DataFrame.from_dict(Descriptores2D)
    Descriptores2D_Dataframe = Descriptores2D_Dataframe.transpose()

    Descriptores2D_Dataframe_subset1 = Descriptores2D_Dataframe.astype('string')
    Descriptores2D_Dataframe_subset1['MQNs_'] = Descriptores2D_Dataframe_subset1['MQNs_'].str.replace("[", "").str.replace(
        "]", "")
    Descriptores2D_Dataframe_subset1 = Descriptores2D_Dataframe_subset1['MQNs_'].str.split(pat=",", expand=True).add_prefix(
        'MQNs_value_').fillna(0)
    Descriptores2D_Dataframe = Descriptores2D_Dataframe.join(Descriptores2D_Dataframe_subset1)
    Descriptores2D_Dataframe = Descriptores2D_Dataframe.drop('MQNs_', axis=1)

    Descriptores2D_Dataframe_subset2 = Descriptores2D_Dataframe.astype('string')
    Descriptores2D_Dataframe_subset2['CalcAUTOCORR2D'] = Descriptores2D_Dataframe_subset2['CalcAUTOCORR2D'].str.replace("[",
                                                                                                                        "").str.replace(
        "]", "")
    Descriptores2D_Dataframe_subset2 = Descriptores2D_Dataframe_subset2['CalcAUTOCORR2D'].str.split(pat=",",
                                                                                                    expand=True).add_prefix(
        'CalcAUTOCORR2D_value_').fillna(0)
    Descriptores2D_Dataframe = Descriptores2D_Dataframe.join(Descriptores2D_Dataframe_subset2)
    Descriptores2D_Dataframe = Descriptores2D_Dataframe.drop('CalcAUTOCORR2D', axis=1)

    Descriptores2D_Dataframe.to_csv("./DescriptorsCSV/Descriptores2D_Dataframe.csv")


    # Exportando base de datos de smiles a una lista

    smiles = pd.read_csv('./SMILES.csv', delimiter=',')
    smiles = smiles.transpose()
    smiles = smiles.values.tolist()

    # Recorriendo la lista para calcular descriptores

    DescriptoresLipinski = {}
    for x in smiles:
        for molecula in x:
            Lipinski_Diccionario = {}
            m = Chem.MolFromSmiles(molecula)
            DescriptoresLipinski_Diccionario = {}

            DescriptoresLipinski_Diccionario['FractionCSP3'] = Lipinski.FractionCSP3(m)
            DescriptoresLipinski_Diccionario['HeavyAtomCount'] = Lipinski.HeavyAtomCount(m)
            DescriptoresLipinski_Diccionario['NHOHCount'] = Lipinski.NHOHCount(m)
            DescriptoresLipinski_Diccionario['NOCount'] = Lipinski.NOCount(m)
            DescriptoresLipinski_Diccionario['NumAliphaticCarbocycles'] = Lipinski.NumAliphaticCarbocycles(m)
            DescriptoresLipinski_Diccionario['NumAliphaticHeterocycles'] = Lipinski.NumAliphaticHeterocycles(m)
            DescriptoresLipinski_Diccionario['NumAromaticCarbocycles'] = Lipinski.NumAromaticCarbocycles(m)
            DescriptoresLipinski_Diccionario['NumHAcceptors'] = Lipinski.NumHAcceptors(m)
            DescriptoresLipinski_Diccionario['NumHDonors'] = Lipinski.NumHDonors(m)
            DescriptoresLipinski_Diccionario['NumHeteroatoms'] = Lipinski.NumHeteroatoms(m)
            DescriptoresLipinski_Diccionario['NumRotatableBonds'] = Lipinski.NumRotatableBonds(m)
            DescriptoresLipinski_Diccionario['NumSaturatedCarbocycles'] = Lipinski.NumSaturatedCarbocycles(m)
            DescriptoresLipinski_Diccionario['NumSaturatedHeterocycles'] = Lipinski.NumSaturatedHeterocycles(m)
            DescriptoresLipinski_Diccionario['NumSaturatedRings'] = Lipinski.NumSaturatedRings(m)
            DescriptoresLipinski_Diccionario['RingCount'] = Lipinski.RingCount(m)

            DescriptoresLipinski[molecula] = DescriptoresLipinski_Diccionario

    # Los descriptores almacenados en un diccionario se convierten a dataframes, se guardan como csv
    DescriptoresLipinski_Dataframe = pd.DataFrame.from_dict(DescriptoresLipinski)
    DescriptoresLipinski_Dataframe = DescriptoresLipinski_Dataframe.transpose()

    DescriptoresLipinski_Dataframe.to_csv("./DescriptorsCSV/DescriptoresLipinski_Dataframe.csv")

if __name__ == "__main__":
    smiles_de_prueba = input("Introduce tu smiles:")
    descriptor3d_de_prueba = calculating_descriptors3d(smiles_de_prueba)
    descriptor3d_de_prueba_transformado = descriptor3d_data_transformation(descriptor3d_de_prueba)
    print(descriptor3d_de_prueba)
    print(descriptor3d_de_prueba_transformado)
