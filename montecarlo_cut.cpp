//
//  Eseguiamo il taglio degli impulsi di interesse dalle simulazioni Montecarlo
//

void montecarlo_cut(){
    // leggiamo il file ed associamolo ad un tree
    TFile *file = new TFile("Jpsi_MC.root");
    TTree *tree = (TTree*) file->Get("data");
    // imponiamo i limiti di massa della JPsi
    const double_t JPsiMassMin = 2.6;
    const double_t JPsiMassMax = 3.4;
    // definiamo gli estremi dell'asse dell'impulso dei muoni nei plot che tracceremo
    double_t PtMin = 0;
    double_t PtMax = 30;
    double_t Nbin = 500;
    // definiamo l'ampiezza dei bin
    double DeltaPt = (PtMax - PtMin) / (1. * Nbin);
    
    // definiamo gli istogrammi per contenere i dati delle quantità da plottare, MC-misurati ed MC-generati
    TH1F *histo1 = new TH1F("mu1_pt", "mu1_pt", Nbin, PtMin, PtMax);
    TH1F *histog1 = new TH1F("mu1_pt_gen", "mu1_pt_gen", Nbin, PtMin, PtMax);
    TH1F *histo2 = new TH1F("mu2_pt", "mu2_pt", Nbin, PtMin, PtMax);
    TH1F *histog2 = new TH1F("mu1_p2_gen", "mu2_pt_gen", Nbin, PtMin, PtMax);
    
    // creiamo vettori di Lorentz per i dati di interesse
    TLorentzVector *mu1 = new TLorentzVector();
    TLorentzVector *gmu1 = new TLorentzVector();
    TLorentzVector *mu2 = new TLorentzVector();
    TLorentzVector *gmu2 = new TLorentzVector();
    TLorentzVector *JPsi = new TLorentzVector();
    TLorentzVector *gJPsi = new TLorentzVector();
    
    // associamo ai vettori di Lorentz i relativi dati, includendo il valore del trigger, la carica della JPsi e la compatibilità montecarlo tra la JPsi e i dati generati
    int JPsiType;
    int MCType, trig, Jq;
    tree->SetBranchAddress("muPosP", &mu1);
    tree->SetBranchAddress("muPosP_Gen", &gmu1);
    tree->SetBranchAddress("muNegP", &mu2);
    tree->SetBranchAddress("muNegP_Gen", &gmu2);
    tree->SetBranchAddress("JpsiP", &JPsi);
    tree->SetBranchAddress("JpsiP_Gen", &gJPsi);
    tree->SetBranchAddress("JpsiType", &JPsiType);
    tree->SetBranchAddress("MCType", &MCType);
    tree->SetBranchAddress("HLT_DoubleMu0", &trig);
    tree->SetBranchAddress("JpsiCharge", &Jq);
    
    // prendiamo il numero totale di eventi registrati
    int ntot = tree->GetEntries();
    int n_scelti = 100000;
    // definiamo array in cui, se soddisferanno le condizioni dell'if, veranno inseriti i dati di impulso e pseudorapidità MC-generati delle tre particelle
    double PtMu1[n_scelti], EtaMu1[n_scelti];
    double PtMu2[n_scelti], EtaMu2[n_scelti];
    double PtJPsi[n_scelti], EtaJPsi[n_scelti];
    
    // stampiamo il numero di eventi totali
    cout << "numero eventi" << ntot << endl;
    
    // definiamo dunque un ciclo for in cui, per ogni giro, viene letto un dato di impulso e pseudorapidità MC-misurato e MC-generato delle tre particelle: se le quantità MC-misurate soddisferanno le condizioni richieste di massa (JPsi) ed impulso (mu1, mu2) + le condizioni di trigger, inseriremo l'evento ed il corrispondente generato negli istogrammi relativi mu[i]_pt e mu[i]_pt_gen; in particolare, per tracciare un plot Pt vs Eta generati, inseriremo i detti dati negli array definiti 7 righe sopra
    int n_letti = 0;
    for(int ev = 0; ev < n_scelti; ev++){
        tree->GetEntry(ev);
        // definiamo gli impulsi
        double Pt1 = mu1->Pt();
        double Eta1 = mu1->PseudoRapidity();
        
        double gPt1 = gmu1->Pt();
        double gEta1 = gmu1->PseudoRapidity();
        
        double Pt2 = mu2->Pt();
        double Eta2 = mu2->PseudoRapidity();
        
        double gPt2 = gmu2->Pt();
        double gEta2 = gmu2->PseudoRapidity();
        
        double PtJP = JPsi->Pt();
        double EtaJP = JPsi->PseudoRapidity();
        
        double gPtJP = gJPsi->Pt();
        double gEtaJP = gJPsi->PseudoRapidity();
        
        // definiamo quindi l'if che ci definisce le condizioni che i successivi eventi devono soddisfare
        if(JPsi->M() > 3 && JPsi->M() < 3.2 && gJPsi->M() > 3 && gJPsi->M() < 3.2 && Jq == 0 && MCType == 0 && (JPsiType == 0 || JPsiType == 1 || JPsiType == 2)){
            if(Pt1 > PtMin && Pt2 > PtMin){
                PtMu1[n_letti] = gPt1;
                EtaMu1[n_letti] = gEta1;
                
                PtMu2[n_letti] = gPt2;
                EtaMu2[n_letti] = gEta2;
                
                PtJPsi[n_letti] = gPtJP;
                EtaJPsi[n_letti] = gEtaJP;
                
                n_letti++;
                histo1->Fill(Pt1);
                histog1->Fill(gPt1);
                histo2->Fill(Pt2);
                histog2->Fill(gPt2);
            }
        }
    }

    cout << "numero eventi letti" << n_letti << endl;
    
    // tracciamo gli istogrammi per gli impulsi MC-misurati e MC-generati
    TCanvas *c = new TCanvas("c");
    c->Divide(1, 2);
    // parte superiore del grafico: mu1
    c->cd(1);
    histog1->SetLineColor(3);
    histo1->SetLineColor(2);
    histog1->Draw("same");
    histo1->Draw("same");
    // parte inferiore del grafico: mu2
    c->cd(2);
    histog2->SetLineColor(3);
    histo2->SetLineColor(2);
    histog2->Draw("same");
    histo2->Draw("same");
    
    // tracciamo il grafico di Pt vs Eta per le tre particelle
    TCanvas *c1 = new TCanvas("c1");
    TGraph *gr1 = new TGraph(n_letti, PtMu1, EtaMu1);
    gr1->SetTitle("Grafico (pT, eta) muone positivo");
    gr1->GetXaxis()->SetTitle("Pt [Gev/c]");
    gr1->GetYaxis()->SetTitle("Pseudorapidity");
    gr1->Draw("AP");
    gr1->SetMarkerStyle(7);
    c1->Print("PtvsEtaPos.root");
    c1->Print("PtvsEtaPos.png");

    TCanvas *c2 = new TCanvas("c2");
    TGraph *gr2 = new TGraph(n_letti, PtMu2, EtaMu2);
    gr2->SetTitle("Grafico (pT, eta) muone negativo");
    gr2->GetXaxis()->SetTitle("Pt [Gev/c]");
    gr2->GetYaxis()->SetTitle("Pseudorapidity");
    gr2->Draw("AP");
    gr2->SetMarkerStyle(7);
    c2->Print("PtvsEtaNeg.root");
    c2->Print("PtvsEtaNeg.png");

    TCanvas *c3 = new TCanvas("c3");
    TGraph *gr3 = new TGraph(n_letti, PtJPsi, EtaJPsi);
    gr3->SetTitle("Grafico (pT, eta) JPsi");
    gr3->GetXaxis()->SetTitle("Pt [Gev/c]");
    gr3->GetYaxis()->SetTitle("Pseudorapidity");
    gr3->Draw("AP");
    gr3->SetMarkerStyle(7);
    c3->Print("PtvsEtaJPsi.root");
    c3->Print("PtvsEtaJPsi.png");
}
