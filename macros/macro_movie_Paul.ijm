// Demander le nom du dossier de sortie
outputDir = getDirectory("Selectionnez ou creez un dossier pour enregistrer les resultats :");

// Vérifier que l'utilisateur a bien sélectionné un dossier
if (outputDir != "") {
    // Demander le nombre de dossiers à traiter
    numFolders = getNumber("Combien de dossiers souhaitez-vous traiter ?", 1);

    if (numFolders > 0) {
        folders = newArray(numFolders);

        // Sélectionner les dossiers à analyser
        for (i = 0; i < numFolders; i++) {
            folders[i] = getDirectory("Selectionnez le dossier " + (i + 1) + " :");
        }

        // Boucle sur chaque dossier sélectionné
        for (k = 0; k < numFolders; k++) {
            dir = folders[k];

            // Découper le chemin en sous-dossiers
            pathParts = split(dir, File.separator);

            // Vérifier qu'il y a au moins deux dossiers dans le chemin
            if (pathParts.length > 2) {
                folderName = pathParts[pathParts.length - 2]; // Avant-dernier dossier
            } else {
                folderName = pathParts[pathParts.length - 1]; // Dernier dossier si structure plus simple
            }
            print("Traitement du dossier : " + dir);

            // Ouvrir la séquence d'images
            run("Image Sequence...", "open=" + dir + " sort");

			run("8-bit");

 			run("Auto Local Threshold", "method=Phansalkar radius=15 parameter_1=0 parameter_2=0 stack");

			run("Set Measurements...", "area center stack redirect=None decimal=3");

			run("Shape Filter", "area=150-1500 area_convex_hull=0-Infinity perimeter=0-Infinity perimeter_convex_hull=0-Infinity feret_diameter=0-Infinity min._feret_diameter=0-Infinity long_side_min._bounding_rect.=0-Infinity short_side_min._bounding_rect.=0-Infinity aspect_ratio=1-Infinity area_to_perimeter_ratio=0-Infinity circularity=0-Infinity elongation=0-1 convexity=0-1 solidity=0-1 num._of_holes=0-1 thinnes_ratio=0-1 contour_temperatur=0-1 fractal_box_dimension=0-2 option->box-sizes=2,3,4,6,8,12,16,32,64 add_to_manager draw_holes black_background fill_results_table stack");

            // Définir le chemin du fichier de sortie avec le nom du dossier source
            resultFilePath = outputDir + folderName + "_Results.csv";

            // Sauvegarde des résultats dans le dossier de sortie
            saveAs("Results", resultFilePath);

            print("Résultats enregistrés sous : " + resultFilePath);

            // Fermeture des fenêtres ouvertes
            selectWindow("Log");
            run("Close");
            selectWindow("Results");
            run("Close");
            selectWindow("ROI Manager");
            run("Close");
        }
    }
}
