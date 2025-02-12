//
//  Tracciamo i plot della popolazione degli impulsi vs rapidità
//

void population(){
    // definiamo i limiti sugli impulsi, dedotti dai grafici del montecarlo di Pt vs Eta, per i quali ha senso considerare gli eventi acquisiti
    double PtMin = 0.;
    double PtMax = 16.;
    double PtMaxJP = 30.;
    // definiamo il numero di bin che vogliamo tracciare nel grafico delle popolazioni
    int Nbin = 35;
    double DeltaPt = (PtMax - PtMin) / Nbin;
    
    // definiamo gli istogrammi per le popolazioni degli eventi per determinati set di valori di Eta e Pt
    TH2D *hgJP = new TH2D("hgJP", "Gen JPsi Eta vs Pt", 50, 0., 2.4, Nbin, PtMin, PtMaxJP);
    TH2D *hg1 = new TH2D("hg1", "Gen Mu+ Eta vs Pt", 50, 0., 2.4, Nbin, PtMin, PtMax);
    TH2D *hg2 = new TH2D("hg2", "Gen Mu- Eta vs Pt", 50, 0., 2.4, Nbin, PtMin, PtMax);
    
    // apriamo il file MC ed accediamo ai dati
    TFile *file = new TFile("Jpsi_MC.root");
    TTree *tree = (TTree*) file->Get("data");
    
    // creiamo i vettori di Lorentz per i dati MC-generati letti dal file
    TLorentzVector *gmu1 = new TLorentzVector;
    TLorentzVector *gmu2 = new TLorentzVector;
    TLorentzVector *gJPsi = new TLorentzVector;
    // creiamo i vettori di Lorentz per i dati MC-misurati letti dal file
    TLorentzVector *mmu1 = new TLorentzVector;
    TLorentzVector *mmu2 = new TLorentzVector;
    TLorentzVector *mJPsi = new TLorentzVector;
    
    int JPsiType;
    int MCType;
    int trig, Jq;
    // indirizziamo i set di dati giusti ai vettori di Lorentz definiti per i dati MC-generati
    tree->SetBranchAddress("muPosP_Gen", &gmu1);
    tree->SetBranchAddress("muNegP_Gen", &gmu2);
    tree->SetBranchAddress("JpsiP_Gen", &gJPsi);
    // e per i dati MC-misurati
    tree->SetBranchAddress("muPosP", &mmu1);
    tree->SetBranchAddress("muNegP", &mmu2);
    tree->SetBranchAddress("JpsiP", &mJPsi);
    // indirizziamo i set di dati giusti alle variabili definite per il tipo di JPsi, per la sua carica e per i trigger
    tree->SetBranchAddress("JpsiType", &JPsiType);
    tree->SetBranchAddress("MCType", &MCType);
    tree->SetBranchAddress("HLT_DoubleMu0", &trig);
    tree->SetBranchAddress("JpsiCharge", &Jq);
    // prendiamo il numero totale di eventi nel tree
    int n_tot = tree->GetEntries();
    // definiamo il ciclo for per riempire i dati che soddisfino le condizioni richieste
    for(int ev = 0; ev < n_tot; ev++){
        // leggiamo l'evento ev-esimo
        tree->GetEntry(ev);
        // prendiamo i dati di pseudorapidità e impulso trasverso MC-generati
        double Eta1 = gmu1->PseudoRapidity();
        double Eta2 = gmu2->PseudoRapidity();
        double EtaJP = fabs(gJPsi->Rapidity());
        double Pt1 = gmu1->Pt();
        double Pt2 = gmu2->Pt();
        double PtJP = gJPsi->Pt();
        
        // riempiamo gli istogrammi relativi ai dati dei muoni e della JPsi
        hgJP->Fill(EtaJP, PtJP);
        hg1->Fill(Eta1, Pt1);
        hg2->Fill(Eta2, Pt2);
        
        // definiamo le variabili booleane che includeranno le condizioni da soddisfare per le pseudo-rapidità e gli impulsi MC-generati dei muoni
        bool taglio1 = (fabs(Eta1) < 0.8 && Pt1 > 1.7) || ((fabs(Eta1) > 0.8 && fabs(Eta1) < 1.7) && Pt1 > 0.9) || (fabs(Eta1) > 1.7 && Pt1 > 0.4);
        bool taglio2 = (fabs(Eta2) < 0.8 && Pt2 > 1.7) || ((fabs(Eta2) > 0.8 && fabs(Eta2) < 1.7) && Pt2 > 0.9) || (fabs(Eta2) > 1.7 && Pt2 > 0.4);
    }
    // tracciamo dunque gli istogrammi: per la J/Psi in particolare tracciamo il grafico con i dati corretti, con le pseudorapidità e gli impulsi trasversi che soddisfino le condizioni richieste
    TCanvas *c1 = new TCanvas("c1");
    hgJP->SetTitle("Grafico (#eta, p_{T}) della J/#Psi");
    hgJP->GetYaxis()->SetTitle("p_{T} [ Gev/c ]");
    hgJP->GetXaxis()->SetTitle("| #eta |");
    hgJP->Draw("COLZ");
    
    TCanvas *c2 = new TCanvas("c2");
    hg1->SetTitle("Grafico (#eta, p_{T}) del #mu^{+}");
    hg1->GetYaxis()->SetTitle("p_{T} [ Gev/c ]");
    hg1->GetXaxis()->SetTitle("| #eta |");
    hg1->Draw("COLZ");
    
    TCanvas *c3 = new TCanvas("c3");
    hg2->SetTitle("Grafico (#eta, p_{T}) del #mu^{-}");
    hg2->GetYaxis()->SetTitle("p_{T} [ Gev/c ]");
    hg2->GetXaxis()->SetTitle("| #eta |");
    hg2->Draw("COLZ");
}
